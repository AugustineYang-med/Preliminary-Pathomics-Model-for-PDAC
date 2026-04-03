import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.calibration import calibration_curve
import os
from skopt import gp_minimize
from skopt.space import Integer
import random

# Set global random seed for reproducibility
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)
random.seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

os.makedirs("./Result", exist_ok=True)

# -------------------------
# 1. Load training/validation data
# -------------------------

df = pd.read_excel("./Data/Pathomics_features_RawData.xlsx")

# Features and labels
X = df.drop(['Patient.ID', 'DFS_days', 'DFS.status', 'OS.Event', 'OS.Time'], axis=1)
y = df['DFS.status']

# -------------------------
# 2. Split training and validation sets
# -------------------------
X_train, X_val, y_train, y_val = train_test_split(
    X, 
    y, 
    test_size=0.3, 
    random_state=42,
    stratify=y
)

# -------------------------
# 3. Data standardization
# -------------------------
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_val = scaler.transform(X_val)

# -------------------------
# 4. Build Transformer model
# -------------------------

class TransformerModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, n_heads, n_layers, output_dim):
        super(TransformerModel, self).__init__()
        self.embedding = nn.Linear(input_dim, hidden_dim)
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d_model=hidden_dim, nhead=n_heads), 
            num_layers=n_layers
        )
        self.fc = nn.Linear(hidden_dim, output_dim)
    
    def forward(self, x):
        x = self.embedding(x)
        x = x.unsqueeze(0)  # Transformer expects a batch-first input
        x = self.transformer(x)
        x = x.squeeze(0)  # Remove the batch dimension
        x = self.fc(x)
        return x

# -------------------------
# 5. Define training process
# -------------------------

def train_and_evaluate(hidden_dim, n_heads, n_layers, device):
    input_dim = X_train.shape[1]
    output_dim = 1  # Binary classification

    model = TransformerModel(input_dim, hidden_dim, n_heads, n_layers, output_dim).to(device)  # Move model to GPU

    criterion = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

    # Create training dataset and data loader
    train_data = TensorDataset(torch.tensor(X_train, dtype=torch.float32).to(device), torch.tensor(y_train.values, dtype=torch.float32).to(device))
    val_data = TensorDataset(torch.tensor(X_val, dtype=torch.float32).to(device), torch.tensor(y_val.values, dtype=torch.float32).to(device))

    train_loader = DataLoader(train_data, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_data, batch_size=64, shuffle=False)

    # Train model
    model.train()
    for epoch in range(10):  # Adjust number of training epochs here
        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs.squeeze(), labels)
            loss.backward()
            optimizer.step()

    # Validate model
    model.eval()
    y_pred_val = []
    with torch.no_grad():
        for inputs, labels in val_loader:
            outputs = model(inputs)
            preds = torch.sigmoid(outputs.squeeze()).round()  # Convert logits to predicted labels
            y_pred_val.extend(preds.cpu().numpy())

    accuracy_val = accuracy_score(y_val, y_pred_val)
    return -accuracy_val  # Bayesian optimization minimizes, so return negative value

def train_final_model(hidden_dim, n_heads, n_layers, device):
    input_dim = X_train.shape[1]
    output_dim = 1
    model = TransformerModel(input_dim, hidden_dim, n_heads, n_layers, output_dim).to(device)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    train_data = TensorDataset(torch.tensor(X_train, dtype=torch.float32).to(device), torch.tensor(y_train.values, dtype=torch.float32).to(device))
    train_loader = DataLoader(train_data, batch_size=64, shuffle=True)
    model.train()
    for epoch in range(10):
        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs.squeeze(), labels)
            loss.backward()
            optimizer.step()
    return model

# -------------------------
# 6. Bayesian optimization hyperparameter search
# -------------------------

# Hyperparameter search space
search_space = [
    Integer(32, 128, name='hidden_dim'),  # Hidden layer dimension
    Integer(2, 8, name='n_heads'),        # Number of attention heads
    Integer(1, 4, name='n_layers')        # Number of Transformer layers
]

# Bayesian optimization: objective is validation accuracy
def objective(params):
    hidden_dim, n_heads, n_layers = params
    # Ensure hidden_dim is divisible by n_heads
    if hidden_dim % n_heads != 0:
        return 1.0  # Return a poor but finite score
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return train_and_evaluate(hidden_dim, n_heads, n_layers, device)

# Run Bayesian optimization
result = gp_minimize(objective, search_space, n_calls=50, random_state=42, verbose=True)

# Print best hyperparameters and results
print("Best parameters:", result.x)
print("Best score:", -result.fun)  # Negate since we returned negative values

# -------------------------
# 7. Train final model with best hyperparameters
# -------------------------
best_hidden_dim, best_n_heads, best_n_layers = result.x

# Initialize model with best hyperparameters
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")  # Check GPU availability
final_model = train_final_model(best_hidden_dim, best_n_heads, best_n_layers, device)

# Train and save final model
torch.save(final_model.state_dict(), 'best_transformer_model_with_optimal_params.pth')

# -------------------------
# 8. Save best model weights
# -------------------------
torch.save(final_model.state_dict(), 'best_transformer_model_with_optimal_params.pth')  # Save best model

# -------------------------
# 9. Evaluate on validation set
# -------------------------
train_data = TensorDataset(torch.tensor(X_train, dtype=torch.float32).to(device), torch.tensor(y_train.values, dtype=torch.float32).to(device))
val_data = TensorDataset(torch.tensor(X_val, dtype=torch.float32).to(device), torch.tensor(y_val.values, dtype=torch.float32).to(device))

train_loader = DataLoader(train_data, batch_size=64, shuffle=True)
val_loader = DataLoader(val_data, batch_size=64, shuffle=False)

final_model.eval()
y_pred_val = []
y_pred_val_prob = []
with torch.no_grad():
    for inputs, labels in val_loader:
        outputs = final_model(inputs)
        preds = torch.sigmoid(outputs.squeeze()).round()  # Convert logits to predicted labels
        prob = torch.sigmoid(outputs.squeeze()).cpu().numpy()  # Move results back to CPU
        y_pred_val.extend(preds.cpu().numpy())
        y_pred_val_prob.extend(prob)

# Validation set evaluation
accuracy_val = accuracy_score(y_val, y_pred_val)
report_val = classification_report(y_val, y_pred_val, output_dict=True)
precision_val = (report_val['1']['precision'] + report_val['0']['precision']) / 2
recall_val = (report_val['1']['recall'] + report_val['0']['recall']) / 2
f1_score_val = (report_val['1']['f1-score'] + report_val['0']['f1-score']) / 2

print(f"Validation Accuracy: {accuracy_val:.4f}")
print(f"Validation Precision: {precision_val:.4f}")
print(f"Validation Recall: {recall_val:.4f}")
print(f"Validation F1-Score: {f1_score_val:.4f}")

# -------------------------
# 10. Generate DCA, confusion matrix, ROC, and calibration curves
# -------------------------
def calculate_net_benefit(y_true, y_pred_prob, thresholds):
    n = len(y_true)
    net_benefits = []
    
    for threshold in thresholds:
        y_pred_thresh = (y_pred_prob >= threshold).astype(int)
        tp = ((y_pred_thresh == 1) & (y_true == 1)).sum()
        fp = ((y_pred_thresh == 1) & (y_true == 0)).sum()
        
        # Calculate Net Benefit
        net_benefit = (tp - fp * (threshold / (1 - threshold))) / n
        net_benefits.append(net_benefit)
    
    return net_benefits

# -------------------------
# 11. Calculate and plot DCA curve
# -------------------------
thresholds_val = np.linspace(0.01, 0.99, 100)
net_benefits_model_val = calculate_net_benefit(y_val, y_pred_val_prob, thresholds_val)
treat_all_net_benefit_val = [((y_val == 1).sum() / len(y_val)) - (t / (1 - t)) for t in thresholds_val]

plt.figure(figsize=(8, 6))
plt.plot(thresholds_val, net_benefits_model_val, label='Transformer Model', color='blue', linewidth=2)
plt.plot(thresholds_val, treat_all_net_benefit_val, label='Treat All', linestyle='-', color='gray')
plt.axhline(y=0, color='black', linestyle='-', label='None')
plt.ylim(-0.6, 0.6)
plt.xlabel('High Risk Threshold')
plt.ylabel('Net Benefit')
plt.title('Decision Curve Analysis (Validation)')
plt.legend(loc='upper right')
plt.grid(True)
plt.savefig('./Result/DCA_Curve_Transformer_val.png', dpi=1000)
plt.close()

# ROC curve
fpr_val, tpr_val, _ = roc_curve(y_val, y_pred_val_prob)
roc_auc_val = roc_auc_score(y_val, y_pred_val_prob)

plt.figure(figsize=(8, 6))
plt.plot(fpr_val, tpr_val, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc_val:.4f})')
plt.plot([0, 1], [0, 1], color='red', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve (Validation)')
plt.legend(loc='lower right')
plt.savefig('./result/ROC_Curve_Transformer_val.png', dpi=1000)
plt.close()

print(f"Val ROC AUC: {roc_auc_val:.4f}")

# Confusion matrix
conf_matrix_val = confusion_matrix(y_val, y_pred_val)
conf_matrix_val_normalized = conf_matrix_val.astype('float') / conf_matrix_val.sum(axis=1)[:, np.newaxis]

plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix_val_normalized, annot=True, fmt='.2f', cmap='Blues',
            xticklabels=['0', '1'], yticklabels=['0', '1'],
            annot_kws={"size": 18})
plt.xlabel('Predicted', fontsize=14)
plt.ylabel('Actual', fontsize=14)
plt.title('Confusion Matrix (Validation)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('./result/Confusion_Matrix_Transformer_val.png', dpi=1000)
plt.close()

# Calibration curve
prob_true_val, prob_pred_val = calibration_curve(y_val, y_pred_val_prob, n_bins=10)

plt.figure(figsize=(8, 6))
plt.plot(prob_pred_val, prob_true_val, marker='o', linewidth=2, label='Transformer')
plt.plot([0, 1], [0, 1], linestyle='--', label='Perfect Calibration')
plt.xlabel('Mean Predicted Probability')
plt.ylabel('True Probability')
plt.title('Calibration Curve (Validation)')
plt.legend(loc='lower right')
plt.grid(True)
plt.savefig('./result/Calibration_Curve_Transformer_val.png', dpi=1000)
plt.close()


# -------------------------
# 12. Calculate and plot test set figures
# -------------------------
# Load test set

df_test = pd.read_excel("./Data/TCGA-PAAD_pathomics_survival_data.xlsx")
X_test = df_test.drop(['sampleid', 'OS.Event', 'OS.Time', 'PFI.Event', 'PFI.Time'], axis=1)
y_test = df_test['PFI.Event']

# Standardize test set data
X_test = scaler.transform(X_test)

# Convert to Tensor
X_test_tensor = torch.tensor(X_test, dtype=torch.float32).to(device)

# Load best model for testing
best_model = TransformerModel(X_train.shape[1], best_hidden_dim, best_n_heads, best_n_layers, 1).to(device)
best_model.load_state_dict(torch.load('best_transformer_model_with_optimal_params.pth'))

best_model.eval()
y_pred_test = []
y_pred_test_prob = []
with torch.no_grad():
    outputs = best_model(X_test_tensor)
    preds = torch.sigmoid(outputs.squeeze()).round()
    prob = torch.sigmoid(outputs.squeeze()).cpu().numpy()
    y_pred_test.extend(preds.cpu().numpy())
    y_pred_test_prob.extend(prob)

accuracy_test = accuracy_score(y_test, y_pred_test)
report_test = classification_report(y_test, y_pred_test, output_dict=True)
precision_test = (report_test['1']['precision'] + report_test['0']['precision']) / 2
recall_test = (report_test['1']['recall'] + report_test['0']['recall']) / 2
f1_score_test = (report_test['1']['f1-score'] + report_test['0']['f1-score']) / 2

print(f"Test Accuracy: {accuracy_test:.4f}")
print(f"Test Precision: {precision_test:.4f}")
print(f"Test Recall: {recall_test:.4f}")
print(f"Test F1-Score: {f1_score_test:.4f}")

# -------------------------
# 13. Generate test set DCA, confusion matrix, ROC, and calibration curves
# -------------------------

# Calculate and plot test set DCA
net_benefits_model_test = calculate_net_benefit(y_test, y_pred_test_prob, thresholds_val)
treat_all_net_benefit_test = [((y_test == 1).sum() / len(y_test)) - (t / (1 - t)) for t in thresholds_val]

plt.figure(figsize=(8, 6))
plt.plot(thresholds_val, net_benefits_model_test, label='Transformer Model', color='blue', linewidth=2)
plt.plot(thresholds_val, treat_all_net_benefit_test, label='Treat All', linestyle='-', color='gray')
plt.axhline(y=0, color='black', linestyle='-', label='None')
plt.ylim(-0.6, 0.6)
plt.xlabel('High Risk Threshold')
plt.ylabel('Net Benefit')
plt.title('Decision Curve Analysis (Test)')
plt.legend(loc='upper right')
plt.grid(True)
plt.savefig('./Result/DCA_Curve_Transformer_test.png', dpi=1000)
plt.close()

# Confusion matrix
conf_matrix_test = confusion_matrix(y_test, y_pred_test)
conf_matrix_test_normalized = conf_matrix_test.astype('float') / conf_matrix_test.sum(axis=1)[:, np.newaxis]

plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix_test_normalized, annot=True, fmt='.2f', cmap='Blues', 
            xticklabels=['0', '1'], yticklabels=['0', '1'],
            annot_kws={"size": 18})
plt.xlabel('Predicted', fontsize=14)
plt.ylabel('Actual', fontsize=14)
plt.title('Confusion Matrix (Test)', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('./Result/Confusion_Matrix_Transformer_test.png', dpi=1000)
plt.close()

# ROC curve
fpr_test, tpr_test, _ = roc_curve(y_test, y_pred_test_prob)
roc_auc_test = roc_auc_score(y_test, y_pred_test_prob)
print(f"Test ROC AUC: {roc_auc_test:.4f}")

plt.figure(figsize=(8, 6))
plt.plot(fpr_test, tpr_test, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc_test:.4f})')
plt.plot([0, 1], [0, 1], color='red', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve (Test)')
plt.legend(loc='lower right')
plt.savefig('./Result/ROC_Curve_Transformer_test.png', dpi=1000)
plt.close()

# Calibration curve
plt.figure(figsize=(8, 6))
prob_true_test, prob_pred_test = calibration_curve(y_test, y_pred_test_prob, n_bins=10)
plt.plot(prob_pred_test, prob_true_test, marker='o', label='Transformer', color='blue')
plt.plot([0, 1], [0, 1], linestyle='--', label='Perfectly Calibrated', color='black')
plt.xlabel('Predicted Probability')
plt.ylabel('True Probability')
plt.title('Calibration Curve (Test)')
plt.legend(loc='upper left')
plt.savefig('./Result/Calibration_Curve_Transformer_test.png', dpi=1000)
plt.close()
