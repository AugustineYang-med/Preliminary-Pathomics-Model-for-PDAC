library(openxlsx)
PKUFH_pathomics_features <- read.xlsx("./Data/PKUFH_pathomics_features_RawData_05-15-2025.xlsx")
PKUFH_pathomics_features <- subset(PKUFH_pathomics_features, select = -c(PatientID, DFS.time, DFS))

PKUFH_features <- PKUFH_pathomics_features[,4:3333]

# 1. Variance filtering
variances <- sapply(PKUFH_features, var, na.rm = TRUE) # Calculate variance for each column
selected_columns <- names(variances[variances > 1.5]) # Select columns with variance > 1.5
filtered_PKUFH_pathomics_features <- PKUFH_features[, selected_columns, drop = FALSE] # Create new data frame with columns having variance > 1.5

# 2. Organize data
PKUFH_survival_OS <- PKUFH_pathomics_features[, 1:3]
PKUFH_pathomics_features <- cbind(PKUFH_survival_OS, filtered_PKUFH_pathomics_features)

# 3. Split dataset
set.seed(123)
library(caret)

PKUFH_pathomics_features$OS <- as.factor(PKUFH_pathomics_features$OS) # Ensure OS is a factor (for 0/1 events)
class(PKUFH_pathomics_features$OS)  # Should be "factor"
table(PKUFH_pathomics_features$OS)  # Check event distribution

# PKUFH dataset stratified sampling (based on OS event)
PKUFH_train_index <- createDataPartition(y = PKUFH_pathomics_features$OS,  # Must be a factor or numeric vector
                                         p = 0.7,                          # Training set ratio (70%)
                                         list = TRUE,                      # Return list format
                                         times = 1                         # Number of samplings (1)
)

# Extract training and test sets
PKUFH_train_data <- PKUFH_pathomics_features[PKUFH_train_index$Resample1, ]  # Training set (70%)
PKUFH_test_data  <- PKUFH_pathomics_features[-PKUFH_train_index$Resample1, ] # Test set (30%)
identical(colnames(PKUFH_train_data), colnames(PKUFH_test_data))

# Check event proportions in training and test sets
prop.table(table(PKUFH_train_data$OS))  # Training set proportion
prop.table(table(PKUFH_test_data$OS))   # Test set proportion

# 4. First perform PCA dimensionality reduction on training set
PKUFH_train_data_features <- subset(PKUFH_train_data, select = -c(ID, OS.time, OS))
PKUFH_train_data_features <- as.data.frame(lapply(PKUFH_train_data_features, function(x) as.numeric(as.character(x)))) # Force convert all columns to numeric type
PKUFH_train_data_features_PCA <- prcomp(PKUFH_train_data_features, scale. = TRUE) # scale. = TRUE standardizes data (mean=0, var=1), recommended for PCA
PKUFH_train_data_features_PCA_100 <- PKUFH_train_data_features_PCA$x[, 1:100, drop = FALSE] # Extract dimensionally reduced data (first 100 principal components)
PKUFH_train_data_features_PCA_100 <- as.data.frame(PKUFH_train_data_features_PCA_100) # Convert results to data frame

summary(PKUFH_train_data_features_PCA)$importance[2, 1:100] # Variance explained proportion
summary(PKUFH_train_data_features_PCA)$importance[3, 1:100] # Cumulative variance explained proportion
PKUFH_train_data_features_PCA$rotation[, 1:100] # Feature contributions to principal components

PKUFH_train_data_features_survival_OS <- subset(PKUFH_train_data, select = c(ID, OS.time, OS))
PKUFH_train_data_features_PCA_100 <- cbind(PKUFH_train_data_features_survival_OS,PKUFH_train_data_features_PCA_100)

# 5. Apply training set PCA transformation to test set
PKUFH_test_data_features <- subset(PKUFH_test_data, select = -c(ID, OS.time, OS))
PKUFH_test_data_features <- as.data.frame(lapply(PKUFH_test_data_features, function(x) as.numeric(as.character(x)))) # Force convert all columns to numeric type
PKUFH_test_data_features_PCA <- predict(PKUFH_train_data_features_PCA, newdata = PKUFH_test_data_features)
PKUFH_test_data_features_PCA_100 <- as.data.frame(PKUFH_test_data_features_PCA) # Convert results to data frame
PKUFH_test_data_features_PCA_100 <- PKUFH_test_data_features_PCA_100[,1:100]
PKUFH_test_data_features_survival_OS <- subset(PKUFH_test_data, select = c(ID, OS.time, OS))
PKUFH_test_data_features_PCA_100 <- cbind(PKUFH_test_data_features_survival_OS,PKUFH_test_data_features_PCA_100)

################################ Model training --------------------------------
############################## Model Construction ------------------------------
# Create a list containing training and test sets
list_train_vali_Data <- list(PKUFH_Trainning_Dataset = PKUFH_train_data_features_PCA_100,
                             PKUFH_Test_Dataset = PKUFH_test_data_features_PCA_100)
features <- colnames(PKUFH_test_data_features_PCA_100)
features <- features[features != c("ID", "OS.time", "OS")]

# View list structure
str(list_train_vali_Data)

# 5. Model construction Important!!!: Model Construction
library(Mime1)

res <- ML_cv5(train_data = list_train_vali_Data$PKUFH_Trainning_Dataset,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = features,
                       mode = 'all',
                       nodesize =5,
                       seed = 123,
                       cores_for_parallel = 60)

cindex_dis_all(res, 
               validate_set = names(list_train_vali_Data)[-1], 
               order = names(list_train_vali_Data), 
               width = 0.35)









col_features <- colnames(PKUFH_pathomics_features_pca)
col_features <- col_features[col_features != "ID"]
col <- c("ID","OS.time", "OS", colnames(PKUFH_pathomics_features_pca))
PKUFH_pathomics_features_pca$ID <- PKUFH_raw_features$ID