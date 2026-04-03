library(openxlsx)
PKUFH_Train_clinical_data <- read.xlsx("./Data/Recurrence/Only Clinical Data/PKUFH_Training_with_Clinical.xlsx")
PKUFH_Test_clinical_data <- read.xlsx("./Data/Recurrence/Only Clinical Data/PKUFH_Test_with_Clinical.xlsx")
TCGA_Test_clinical_data <- read.xlsx("./Data/Recurrence/Only Clinical Data/TCGA_Test_with_Clinical.xlsx")

identical(colnames(PKUFH_Train_clinical_data), colnames(PKUFH_Train_clinical_data))
PKUFH_Clinical_data <- rbind(PKUFH_Train_clinical_data, PKUFH_Test_clinical_data)

PKUFH_Train_group <- subset(PKUFH_train_data, select = c(ID,Var))
PKUFH_Test_group <- subset(PKUFH_test_data, select = c(ID,Var))
TCGA_Test_group <- subset(TCGA_Test_Data, select = c(ID,Var))

PKUFH_Train_dataset_only_clinical <- merge(PKUFH_Train_group, PKUFH_Clinical_data, by.x = "ID", by.y = "ID")
PKUFH_Test_dataset_only_clinical <- merge(PKUFH_Test_group, PKUFH_Clinical_data, by.x = "ID", by.y = "ID")
TCGA_Test_dataset_only_clinical <- merge(TCGA_Test_group, TCGA_Test_clinical_data, by.x = "ID", by.y = "ID")

PKUFH_Train_dataset_only_clinical <- subset(PKUFH_Train_dataset_only_clinical, select = -c(OS.time,OS,DFS.time,DFS))
PKUFH_Test_dataset_only_clinical <- subset(PKUFH_Test_dataset_only_clinical, select = -c(OS.time,OS,DFS.time,DFS))
TCGA_Test_dataset_only_clinical <- subset(TCGA_Test_dataset_only_clinical, select = -c(OS.time,OS,DFS.time,DFS))

cols <- colnames(PKUFH_Train_dataset_only_clinical)
TCGA_Test_dataset_only_clinical <- TCGA_Test_dataset_only_clinical[,cols]
identical(colnames(PKUFH_Train_dataset_only_clinical), colnames(PKUFH_Test_dataset_only_clinical))
identical(colnames(PKUFH_Test_dataset_only_clinical), colnames(TCGA_Test_dataset_only_clinical))

PKUFH_Train_dataset_only_clinical$Pathology <- ifelse(PKUFH_Train_dataset_only_clinical$Pathology == "G1", 1,
                                                      ifelse(PKUFH_Train_dataset_only_clinical$Pathology == "G2", 2, 3))
PKUFH_Train_dataset_only_clinical$R0.resection <- ifelse(PKUFH_Train_dataset_only_clinical$R0.resection == "Yes", 1,
                                                         ifelse(PKUFH_Train_dataset_only_clinical$R0.resection == "No", 0, NA))
PKUFH_Train_dataset_only_clinical$Tumor.site <- ifelse(PKUFH_Train_dataset_only_clinical$Tumor.site == "Body/Tail", 1,
                                                       ifelse(PKUFH_Train_dataset_only_clinical$Tumor.site == "Head/Neck", 0, NA))
PKUFH_Train_dataset_only_clinical$AJCC_8th <- ifelse(PKUFH_Train_dataset_only_clinical$AJCC_8th == "<IIB", 0,
                                                     ifelse(PKUFH_Train_dataset_only_clinical$AJCC_8th == ">=IIB", 1, NA))

TCGA_Test_dataset_only_clinical$Pathology <- ifelse(TCGA_Test_dataset_only_clinical$Pathology == "G1", 1,
                                                    ifelse(TCGA_Test_dataset_only_clinical$Pathology == "G2", 2, 3))
TCGA_Test_dataset_only_clinical$R0.resection <- ifelse(TCGA_Test_dataset_only_clinical$R0.resection == "Yes", 1,
                                                       ifelse(TCGA_Test_dataset_only_clinical$R0.resection == "No", 0, NA))
TCGA_Test_dataset_only_clinical$Tumor.site <- ifelse(TCGA_Test_dataset_only_clinical$Tumor.site == "Body/Tail", 1,
                                                     ifelse(TCGA_Test_dataset_only_clinical$Tumor.site == "Head/Neck", 0, NA))
TCGA_Test_dataset_only_clinical$AJCC_8th <- ifelse(TCGA_Test_dataset_only_clinical$AJCC_8th == "<IIB", 0,
                                                   ifelse(TCGA_Test_dataset_only_clinical$AJCC_8th == ">=IIB", 1, NA))

PKUFH_Test_dataset_only_clinical$Pathology <- ifelse(PKUFH_Test_dataset_only_clinical$Pathology == "G1", 1,
                                                     ifelse(PKUFH_Test_dataset_only_clinical$Pathology == "G2", 2, 3))
PKUFH_Test_dataset_only_clinical$R0.resection <- ifelse(PKUFH_Test_dataset_only_clinical$R0.resection == "Yes", 1,
                                                        ifelse(PKUFH_Test_dataset_only_clinical$R0.resection == "No", 0, NA))
PKUFH_Test_dataset_only_clinical$Tumor.site <- ifelse(PKUFH_Test_dataset_only_clinical$Tumor.site == "Body/Tail", 1,
                                                      ifelse(PKUFH_Test_dataset_only_clinical$Tumor.site == "Head/Neck", 0, NA))
PKUFH_Test_dataset_only_clinical$AJCC_8th <- ifelse(PKUFH_Test_dataset_only_clinical$AJCC_8th == "<IIB", 0,
                                                    ifelse(PKUFH_Test_dataset_only_clinical$AJCC_8th == ">=IIB", 1, NA))

save(PKUFH_Train_dataset_only_clinical, PKUFH_Test_dataset_only_clinical, TCGA_Test_dataset_only_clinical,
     file = "./Data/Recurrence/Only Clinical Data/Three Datasets with Only Clinical Data_features Encoded.RData")

colnames(PKUFH_Train_dataset_only_clinical)[colnames(PKUFH_Train_dataset_only_clinical) == "AJCC_8th"] <- "AJCC"
colnames(PKUFH_Test_dataset_only_clinical)[colnames(PKUFH_Test_dataset_only_clinical) == "AJCC_8th"] <- "AJCC"
colnames(TCGA_Test_dataset_only_clinical)[colnames(TCGA_Test_dataset_only_clinical) == "AJCC_8th"] <- "AJCC"

PKUFH_Train_dataset_only_clinical$Pathology <- as.factor(PKUFH_Train_dataset_only_clinical$Pathology)
PKUFH_Train_dataset_only_clinical$R0.resection <- as.factor(PKUFH_Train_dataset_only_clinical$R0.resection)
PKUFH_Train_dataset_only_clinical$Tumor.site <- as.factor(PKUFH_Train_dataset_only_clinical$Tumor.site)
PKUFH_Train_dataset_only_clinical$AJCC <- as.factor(PKUFH_Train_dataset_only_clinical$AJCC)
PKUFH_Train_dataset_only_clinical$pT.stage <- as.factor(PKUFH_Train_dataset_only_clinical$pT.stage)
PKUFH_Train_dataset_only_clinical$pN.stage <- as.factor(PKUFH_Train_dataset_only_clinical$pN.stage)

TCGA_Test_dataset_only_clinical$Var <- as.factor(TCGA_Test_dataset_only_clinical$Var)

PKUFH_Test_dataset_only_clinical$Pathology <- as.factor(PKUFH_Test_dataset_only_clinical$Pathology)
PKUFH_Test_dataset_only_clinical$R0.resection <- as.factor(PKUFH_Test_dataset_only_clinical$R0.resection)
PKUFH_Test_dataset_only_clinical$Tumor.site <- as.factor(PKUFH_Test_dataset_only_clinical$Tumor.site)
PKUFH_Test_dataset_only_clinical$AJCC <- as.factor(PKUFH_Test_dataset_only_clinical$AJCC)
PKUFH_Test_dataset_only_clinical$pT.stage <- as.factor(PKUFH_Test_dataset_only_clinical$pT.stage)
PKUFH_Test_dataset_only_clinical$pN.stage <- as.factor(PKUFH_Test_dataset_only_clinical$pN.stage)


TCGA_Test_dataset_only_clinical$Pathology <- as.factor(TCGA_Test_dataset_only_clinical$Pathology)
TCGA_Test_dataset_only_clinical$R0.resection <- as.factor(TCGA_Test_dataset_only_clinical$R0.resection)
TCGA_Test_dataset_only_clinical$Tumor.site <- as.factor(TCGA_Test_dataset_only_clinical$Tumor.site)
TCGA_Test_dataset_only_clinical$AJCC <- as.factor(TCGA_Test_dataset_only_clinical$AJCC)
TCGA_Test_dataset_only_clinical$pT.stage <- as.factor(TCGA_Test_dataset_only_clinical$pT.stage)
TCGA_Test_dataset_only_clinical$pN.stage <- as.factor(TCGA_Test_dataset_only_clinical$pN.stage)

PKUFH_train_data <- subset(PKUFH_train_data, select = -c(Var))
PKUFH_Train_pathomics_clinical <- merge(PKUFH_Train_dataset_only_clinical, PKUFH_train_data, by.x = "ID", by.y = "ID")

PKUFH_test_data <- subset(PKUFH_test_data, select = -c(Var))
PKUFH_Test_pathomics_clinical <- merge(PKUFH_Test_dataset_only_clinical, PKUFH_test_data, by.x = "ID", by.y = "ID")

TCGA_Test_Data <- subset(TCGA_Test_Data, select = -c(Var))
TCGA_Test_pathomics_clinical <- merge(TCGA_Test_dataset_only_clinical, TCGA_Test_Data, by.x = "ID", by.y = "ID")


################################ Model Training --------------------------------
############################## Model Construction ------------------------------
# Create a list containing training and test sets
list_train_vali_Data <- list(PKUFH_Training_Dataset = PKUFH_Train_pathomics_clinical,
                             PKUFH_Test_Dataset = PKUFH_Test_pathomics_clinical,
                             TCGA_Test_Dataset = TCGA_Test_pathomics_clinical)
# View list structure
str(list_train_vali_Data)

col <- colnames(TCGA_Test_pathomics_clinical)
col_features <- col[col != c("ID","Var")]
# 5. Important!!!: Model Construction
library(Mime1)

res_recurrence_with_clinical_data <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$PKUFH_Training_Dataset,
                                                              list_train_vali_Data = list_train_vali_Data,
                                                              candidate_genes = col_features,
                                                              methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                                              seed = 123,
                                                              cores_for_parallel = 20)
auc_vis_category_all(res_recurrence_with_clinical_data,
                     dataset = c("PKUFH_Training_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"),
                     order= c("PKUFH_Training_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"))

pdf("./Result/Plot/AUC_Recurrence_with_clinical_data_c.pdf", width = 20, height = 16)  # Adjust width and height as needed
plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res_recurrence_with_clinical_data,model_name = i,dataset = c("PKUFH_Training_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"),
                                   order= c("PKUFH_Training_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)
dev.off()

save(res_recurrence_with_clinical_data, file = "./Result/RData/Recurrence_Model_with_clinicaldata.RData")






