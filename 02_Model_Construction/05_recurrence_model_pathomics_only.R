library(openxlsx)
PKUFH_recurrence_tag <- read.xlsx("./Data/Recurrence/Recurrence OR No Recurrence_tag.xlsx")
TCGA_Test_Data <- read.xlsx("./Result/PCA/PCA100/TCGA_data_PFI_Model.xlsx")
TCGA_Test_Data <- subset(TCGA_Test_Data, select = -c(PFI.time))

colnames(TCGA_Test_Data)[colnames(TCGA_Test_Data) == "PFI"] <- "Var"

PKUFH_PCA_results <- subset(PKUFH_PCA_results,select = -c(OS.time,OS))
col <- colnames(PKUFH_PCA_results)
col <- col[col != "ID"]
col <- c("ID","Var",col)

PKUFH_Recurrence_PCA100_Results <- merge(recurrence_tag, PKUFH_PCA_results, by.x = "ID", by.y = "ID")

library(dplyr)
PKUFH_Recurrence_PCA100_Results <- PKUFH_Recurrence_PCA100_Results %>% 
  mutate(Var = recode(Var, "Yes" = "Y", "No" = "N"))

TCGA_Test_Data <- TCGA_Test_Data %>% 
  mutate(Var = recode(Var, "1" = "Y", "0" = "N"))



# 2. Split dataset
set.seed(123)
library(caret)

PKUFH_Recurrence_PCA100_Results$Var <- as.factor(PKUFH_Recurrence_PCA100_Results$Var) # Ensure Var is a factor (for 0/1 events)
class(PKUFH_Recurrence_PCA100_Results$Var)  # Should be "factor"
table(PKUFH_Recurrence_PCA100_Results$Var)  # Check event distribution

# PKUFH dataset stratified sampling (based on OS event)
PKUFH_train_index <- createDataPartition(y = PKUFH_Recurrence_PCA100_Results$Var,  # Must be a factor or numeric vector
                                         p = 0.7,                          # Training set ratio (70%)
                                         list = TRUE,                      # Return list format
                                         times = 1                         # Number of samplings (1)
)

# Extract training and test sets
PKUFH_train_data <- PKUFH_Recurrence_PCA100_Results[PKUFH_train_index$Resample1, ]  # Training set (70%)
PKUFH_test_data  <- PKUFH_Recurrence_PCA100_Results[-PKUFH_train_index$Resample1, ] # Test set (30%)
identical(colnames(PKUFH_train_data), colnames(PKUFH_test_data))

# Check event proportions in training and test sets
prop.table(table(PKUFH_train_data$Var))  # Training set proportion
prop.table(table(PKUFH_test_data$Var))   # Test set proportion
save(PKUFH_Recurrence_PCA100_Results, PKUFH_train_data, PKUFH_test_data, file = "./Data/PKUFH_train_validation_Recurrence Data.RData")

################################ Model Training --------------------------------
############################## Model Construction ------------------------------
# Create a list containing training and test sets
list_train_vali_Data <- list(PKUFH_Training_Dataset = PKUFH_train_data,
                             PKUFH_Test_Dataset = PKUFH_test_data,
                             TCGA_Test_Dataset = TCGA_Test_Data)
# View list structure
str(list_train_vali_Data)

col <- colnames(TCGA_Test_Data)
col_features <- col[col != c("ID","Var")]
# 5. Important!!!: Model Construction
library(Mime1)

res.recurrence <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$PKUFH_Trainning_Dataset,
                                           list_train_vali_Data = list_train_vali_Data,
                                           candidate_genes = col_features,
                                           methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                           seed = 123,
                                           cores_for_parallel = 20)
auc_vis_category_all(res.recurrence,
                     dataset = c("PKUFH_Trainning_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"),
                     order= c("PKUFH_Trainning_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"))

pdf("./Result/Plot/AUC_Recurrence_without_clinical_data_combined.pdf", width = 20, height = 16)  # Adjust width and height as needed
plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.recurrence,model_name = i,dataset = c("PKUFH_Trainning_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"),
                                   order= c("PKUFH_Trainning_Dataset","PKUFH_Test_Dataset", "TCGA_Test_Dataset"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)
dev.off()

save(res.recurrence, file = "./Result/RData/Recurrence_Model_without_clinicaldata.RData")






