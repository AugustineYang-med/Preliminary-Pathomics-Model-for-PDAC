# PCPS_DL: A Deep Learning Pathomic Signature Predicts Survival and Recurrence in Pancreatic Ductal Adenocarcinoma

This repository contains the analysis code for the paper:

> **A Deep Learning Pathomic Signature Predicts Survival and Recurrence in Pancreatic Ductal Adenocarcinoma Independent of CA19-9**
>
> Bohan Yang, Xinxin Liu, Jianlei Wei, Lizhi Xu, Jixin Zhang, Minghe Fan, Yukun Hou, Fusheng Zhang, Yiran Chen, Weikang Liu, Yangyang Li, Zhuo Liu, Kai Chen\*, Xiaodong Tian\*, Yinmo Yang\*, Yongsu Ma\*
>
> Department of Hepatobiliary and Pancreatic Surgery, Peking University First Hospital

## Overview

We developed and validated **PCPS_DL** (PDAC CellProfiler Pathomic Signature - Deep Learning), a deep learning-based pathomic prognostic model for pancreatic ductal adenocarcinoma (PDAC) using routinely acquired H&E-stained whole-slide images (WSIs).

### Key Findings

- **PCPS_DL** robustly predicted overall survival (OS) and recurrence-free survival (RFS) across internal and external cohorts (C-indices: 0.64 and 0.61 internally, 0.62 and 0.61 externally)
- Effectively stratified **CA19-9-negative** patients into high- and low-risk groups (OS: HR = 9.26, P < 0.001; RFS: HR = 5.02, P < 0.001)
- High-risk scores were associated with an immunosuppressive microenvironment and basal-like epithelial cell state
- The model provides biologically interpretable prognostic biomarkers beyond conventional serological markers

### Study Design

```
PKUFH Cohort (n=252)                    TCGA Cohort (n=183)
  ├── Training (n=177)                    └── External Validation
  └── Internal Test (n=75)

H&E WSIs → Tiling → HALO Segmentation → CellProfiler Feature Extraction
       → PCA Dimensionality Reduction → DeepSurv Survival Model
       → Risk Stratification → Biological Interpretation
```

## Repository Structure

```
.
├── 01_WSI_Preprocessing/           # Whole-slide image tiling
│   ├── split_tif.py                # Tile TIF format WSIs with tissue filtering
│   └── split_svs.py                # Tile SVS format WSIs using OpenSlide
│
├── 02_Model_Construction/          # Prognostic model building
│   ├── 01_OS_model_construction.R         # OS prediction model (train/test split)
│   ├── 02_DFS_model_pathomics_only.R      # DFS model using pathomics features only
│   ├── 03_DFS_model_with_clinical.R       # DFS model with clinical + pathomics
│   ├── 04_recurrence_data_preparation.R   # Recurrence model data preprocessing
│   ├── 05_recurrence_model_pathomics_only.R  # Recurrence model (pathomics only)
│   ├── 06_recurrence_model_with_clinical.R   # Recurrence model (pathomics + clinical)
│   ├── 07_recurrence_model_clinical_only.R   # Recurrence model (clinical only)
│   ├── Deep_Learning/
│   │   ├── DeepSurv_model.py              # DeepSurv (CoxPH) survival model with
│   │   │                                  # Bayesian hyperparameter optimization
│   │   └── Transformer_model.py           # Transformer-based classification model
│   └── ML_Framework/
│       ├── ML_ensemble_cv5.R              # Mime ML ensemble (5-fold CV)
│       └── ML_ensemble_cv10.R             # Mime ML ensemble (10-fold CV)
│
├── 03_Differential_Expression/     # DEG analysis between risk groups
│   ├── 01_DESeq2_analysis.R               # DEG identification using DESeq2
│   ├── 02_edgeR_analysis.R                # DEG identification using edgeR
│   ├── 03_limma_analysis.R                # DEG identification using limma
│   └── 07_DEG_visualization.R             # Volcano plots and survival curves for DEGs
│
├── README.md
├── requirements.txt
├── .gitignore
└── LICENSE
```

## Methods

### 1. WSI Preprocessing

H&E-stained slides were scanned at 20x magnification and converted to SVS format. WSIs were tiled into **1024 x 1024 pixel patches** at 0.504 um/pixel resolution. A HALO-based random forest classifier segmented tissue regions, retaining patches with **>=20% tumor content** and **<5% blank space**.

### 2. Feature Extraction

**CellProfiler v4.2.8** was used to extract 3,447 quantitative pathomics features per patch:
- **Nuclear features** (570): morphology, texture, intensity
- **Cytoplasmic features** (569): granularity, intensity distribution
- **Cell membrane features** (570): shape, texture
- **Image-level features** (1,738): quality, colocalization, spatial distribution

After excluding 117 irrelevant features, 3,330 features were retained. Variance filtering (variance >= 1.5) yielded 883 features, which were reduced to **100 principal components** via PCA.

### 3. Model Construction

The **DeepSurv** framework (MLP-based Cox proportional hazards deep neural network) was used to build prognostic models for OS and RFS prediction:
- Activation: ReLU
- Optimizer: Adam
- Hyperparameter tuning: Bayesian search (>40,000 iterations)
- Validation: 5-fold cross-validation
- Comparison: 15+ ML algorithms via the Mime ensemble framework

### 4. Biological Interpretation

- **SHAP analysis**: Feature importance and contribution visualization
- **Differential expression**: DESeq2 + edgeR for DEG identification between risk groups
- **GO/KEGG enrichment**: Functional pathway analysis of DEGs
- **CIBERSORTx**: Immune cell deconvolution from bulk RNA-seq
- **EcoTyper**: Carcinoma ecotype identification and cell state analysis
- **mIHC validation**: Multiplex immunohistochemistry for lymphoid and myeloid panels

## Requirements

### Python (>= 3.9)

```
torch >= 1.12
pycox >= 0.3.0
scikit-learn >= 1.0
scikit-optimize >= 0.9
numpy
pandas
matplotlib
Pillow
openslide-python
histolab >= 0.6.0
```

### R (>= 4.2.0)

```
caret, glmnet, randomForestSRC, gbm, CoxBoost, survivalsvm
survival, survminer, timeROC, rms
DESeq2, edgeR, limma, clusterProfiler
readxl, openxlsx, dplyr, ggplot2
```

See `requirements.txt` for the complete list.

## Data Availability

- **TCGA-PAAD data**: Available from [UCSC Xena](https://xena.ucsc.edu/) and [GDC Data Portal](https://portal.gdc.cancer.gov/)
- **PKUFH cohort data**: Due to patient privacy, access to de-identified data for non-commercial research can be requested from the corresponding author (B.H.Y., ybh000814@163.com). All reasonable requests will be addressed within 14 working days.
- **CIBERSORTx**: https://cibersortx.stanford.edu/
- **EcoTyper**: https://github.com/digitalcytometry/ecotyper

## Software & Tools

| Tool | Version | Purpose |
|------|---------|---------|
| CellProfiler | 4.2.8 | Pathomics feature extraction |
| HALO | 4.0.5107.407 | Tissue segmentation & mIHC quantification |
| DeepSurv (pycox) | 0.3.0 | Deep learning survival model |
| Mime | 0.0.0.9 | ML ensemble framework |
| DESeq2 | 1.38.3 | Differential expression analysis |
| edgeR | 3.36.0 | Differential expression analysis |
| clusterProfiler | 4.10.1 | GO/KEGG enrichment analysis |
| GraphPad Prism | 10.1.2 | Statistical analysis & visualization |
| SHAP | 0.48.0 | Model interpretability |

## Citation

If you use this code, please cite:

```
Yang B, Liu X, Wei J, Xu L, Zhang J, et al. A Deep Learning Pathomic Signature
Predicts Survival and Recurrence in Pancreatic Ductal Adenocarcinoma Independent
of CA19-9. Cancer Letters. 2025.
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Bohan Yang** - ybh000814@163.com
- **Yongsu Ma** - mys870311@sina.com
- **Yinmo Yang** - yangyinmosci@bjmu.edu.cn
- **Xiaodong Tian** - tianxiaodong@pkufh.com
- **Kai Chen** - drchenkai@pku.edu.cn

Department of Hepatobiliary and Pancreatic Surgery, Peking University First Hospital, Beijing, China
