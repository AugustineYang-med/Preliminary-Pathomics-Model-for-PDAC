import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# For preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn_pandas import DataFrameMapper
import torch # For building the networks 
import torchtuples as tt # Some useful functions
# from pycox.datasets import metabric
from pycox.models import LogisticHazard
# from pycox.models import PMF
# from pycox.models import DeepHitSingle
from pycox.evaluation import EvalSurv

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn_pandas import DataFrameMapper

import torch
import torchtuples as tt

from pycox.datasets import metabric
from pycox.models import CoxPH
from pycox.evaluation import EvalSurv

from lifelines import KaplanMeierFitter

import os

import optuna
import json

# We also set some seeds to make this reproducable.
# Note that on gpu, there is still some randomness.
np.random.seed(1234)
_ = torch.manual_seed(123)

# Read internal and external test datasets from xlsx files
df_all = pd.read_excel('./Data/Pathomics_features_RawData.xlsx')

# Shuffle the data
np.random.seed(1234)
df_all = df_all.sample(frac=1, random_state=1234).reset_index(drop=True)

n_total = len(df_all)
n_train = int(n_total * 0.6)
n_val = int(n_total * 0.2)
n_test = n_total - n_train - n_val  # Remaining samples for test

df_train = df_all.iloc[:n_train]
df_val = df_all.iloc[n_train:n_train + n_val]
df_test = df_all.iloc[n_train + n_val:]

cols_standardize = df_train.columns[5:]
cols_leave = []
standardize = [([col], StandardScaler()) for col in cols_standardize]
leave = [(col, None) for col in cols_leave]

x_mapper = DataFrameMapper(standardize + leave)

x_train = x_mapper.fit_transform(df_train).astype('float32')
x_val = x_mapper.transform(df_val).astype('float32')
x_test = x_mapper.transform(df_test).astype('float32')

get_target = lambda df: (df['OS.Time'].values, df['OS.Event'].values)
y_train = get_target(df_train)
y_val = get_target(df_val)
durations_test, events_test = get_target(df_test)
val = x_val, y_val

in_features = x_train.shape[1]
out_features = 1
batch_norm = True
dropout = 0.1
output_bias = False

df_internal = pd.concat([df_train, df_val])
x_internal = x_mapper.transform(df_internal).astype('float32')
durations_internal, events_internal = get_target(df_internal)

results = []

# Bayesian optimization objective function
# Save model for each trial, evaluate on validation set only

def objective(trial):
    num_nodes = [trial.suggest_int('num_nodes1', 16, 128), trial.suggest_int('num_nodes2', 16, 128)]
    dropout = trial.suggest_float('dropout', 0.0, 0.5)
    lr = trial.suggest_float('lr', 1e-4, 1e-2, log=True)
    batch_size = trial.suggest_categorical('batch_size', [64, 128, 256])
    
    net = tt.practical.MLPVanilla(in_features, num_nodes, out_features, batch_norm, dropout, output_bias=output_bias)
    model = CoxPH(net, tt.optim.Adam)
    model.optimizer.set_lr(lr)
    
    callbacks = [tt.callbacks.EarlyStopping()]
    verbose = False
    epochs = 512
    val = x_val, y_val
    try:
        log = model.fit(x_train, y_train, batch_size, epochs, callbacks, verbose, val_data=val, val_batch_size=batch_size)
        model.compute_baseline_hazards(x_train, y_train)
        
        # Training set C-index
        surv_train = model.predict_surv_df(x_train)
        ev_train = EvalSurv(surv_train, y_train[0], y_train[1], censor_surv='km')
        c_index_train = ev_train.concordance_td('antolini')
        
        # Validation set C-index
        surv_val = model.predict_surv_df(x_val)
        ev_val = EvalSurv(surv_val, y_val[0], y_val[1], censor_surv='km')
        c_index_val = ev_val.concordance_td('antolini')
        
        # Test set C-index
        surv_test = model.predict_surv_df(x_test)
        ev_test = EvalSurv(surv_test, durations_test, events_test, censor_surv='km')
        c_index_test = ev_test.concordance_td('antolini')
        
        # Save model
        model_file = f'model_trial_{trial.number}.pt'
        model.save_net(model_file)
        # Log results
        results.append({
            'params': trial.params,
            'train_cindex': c_index_train,
            'val_cindex': c_index_val,
            'test_cindex': c_index_test,
            'model_file': model_file
        })
        return c_index_val
    except Exception as e:
        # If training fails, return very low score
        results.append({'params': trial.params, 'val_cindex': 0.0, 'model_file': None, 'error': str(e)})
        return 0.0

# Start Bayesian optimization
study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=100)

print('All results:')
for r in results:
    print(r)

# Save as JSON file
with open('results.json', 'w') as f:
    json.dump(results, f, indent=2)

# Select model with highest validation C-index
best_result = max(results, key=lambda x: x['val_cindex'])
best_model_file = best_result['model_file']
print(f'Best model file: {best_model_file}, val_cindex: {best_result["val_cindex"]}')

# Load best model
best_net = tt.practical.MLPVanilla(in_features, [best_result['params']['num_nodes1'], best_result['params']['num_nodes2']], out_features, batch_norm, best_result['params']['dropout'], output_bias=output_bias)
best_model = CoxPH(best_net, tt.optim.Adam)
best_model.load_net(best_model_file)
best_model.compute_baseline_hazards(x_train, y_train)

# Evaluate best model on external test set
surv = best_model.predict_surv_df(x_test)
time_point = surv.index[len(surv.index) // 2]   # Example: Middle time point.  Choose a relevant time.
risk_scores = -np.log(surv.loc[time_point]).values
durations_test, events_test = get_target(df_test)

# Calculate median risk score
median_risk = np.median(risk_scores)

# Group assignment
high_risk_group = risk_scores >= median_risk
low_risk_group = risk_scores < median_risk

# 6. Fit Kaplan-Meier curves for each group
kmf_high = KaplanMeierFitter()
kmf_low = KaplanMeierFitter()

kmf_high.fit(durations_test[high_risk_group], event_observed=events_test[high_risk_group], label='High Risk')
kmf_low.fit(durations_test[low_risk_group], event_observed=events_test[low_risk_group], label='Low Risk')

# 7. Plot survival curves
plt.figure(figsize=(8, 6))
kmf_high.plot_survival_function()
kmf_low.plot_survival_function()
plt.title('Survival Curves by Risk Group (Best Model)')
plt.xlabel('Time')
plt.ylabel('Survival Probability')
plt.grid(True)
plt.savefig('./Result/survival_curves_by_risk_group.png')
plt.close()

from lifelines.statistics import logrank_test

y_high_risk = df_test[high_risk_group]
y_low_risk = df_test[low_risk_group]
# Step 7: Perform log-rank test to compare survival between groups
logrank_results = logrank_test(y_high_risk['OS.Time'], y_low_risk['OS.Time'],
                       event_observed_A=y_high_risk['OS.Event'], event_observed_B=y_low_risk['OS.Event'])
print('Logrank test p-value:', logrank_results.p_value)

surv.iloc[:, :5].plot()
plt.ylabel('S(t | x)')
_ = plt.xlabel('Time')
plt.savefig('./Result/survival_probability_examples.png')
plt.close()

ev = EvalSurv(surv, durations_test, events_test, censor_surv='km')

time_grid = np.linspace(durations_test.min(), durations_test.max(), 100)
_ = ev.brier_score(time_grid).plot()
plt.savefig('./Result/brier_score.png')
plt.close()
