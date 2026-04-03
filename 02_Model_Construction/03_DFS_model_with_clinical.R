set.seed(123)  # Set random seed for reproducibility
library(openxlsx)
PKUFH_Train_Dataset <- read.xlsx("./Result/OS Prognostic Data/data with only clinical data/PKUFH_Training_with_Clinical.xlsx")
PKUFH_Test_Dataset <- read.xlsx("./Result/OS Prognostic Data/data with only clinical data/PKUFH_Test_with_Clinical.xlsx")
TCGA_Test_Dataset <- read.xlsx("./Result/OS Prognostic Data/data with only clinical data/TCGA_Test_with_Clinical.xlsx")

colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS"] <- "OS"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS"] <- "OS"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "PFI.time"] <- "OS.time"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "PFI"] <- "OS"

cols <- colnames(TCGA_Test_Dataset)
cols <- cols[!cols %in% c("ID", "DFS.time", "DFS", "OS.time", "OS")]
cols <- c("ID", "OS.time", "OS", "DFS.time", "DFS", cols)
PKUFH_Train_Dataset <- PKUFH_Train_Dataset[,cols]
PKUFH_Test_Dataset <- PKUFH_Test_Dataset[,cols]
TCGA_Test_Dataset <- TCGA_Test_Dataset[,cols]

features <- colnames(TCGA_Test_Dataset)
features <- features[!features %in% c("ID", "DFS.time", "DFS", "OS.time", "OS")]
features
# [1] "Pathology"    "R0.resection" "pT.stage"     "pN.stage"     "Tumor.site"   "AJCC_8th" 
# Assuming your data frame is named df
PKUFH_Train_Dataset$Pathology <- ifelse(PKUFH_Train_Dataset$Pathology == "G1", 1,
                                 ifelse(PKUFH_Train_Dataset$Pathology == "G2", 2, 3))
PKUFH_Train_Dataset$R0.resection <- ifelse(PKUFH_Train_Dataset$R0.resection == "Yes", 1,
                                    ifelse(PKUFH_Train_Dataset$R0.resection == "No", 0, NA))
PKUFH_Train_Dataset$Tumor.site <- ifelse(PKUFH_Train_Dataset$Tumor.site == "Body/Tail", 1,
                                  ifelse(PKUFH_Train_Dataset$Tumor.site == "Head/Neck", 0, NA))
PKUFH_Train_Dataset$AJCC_8th <- ifelse(PKUFH_Train_Dataset$AJCC_8th == "<IIB", 0,
                                ifelse(PKUFH_Train_Dataset$AJCC_8th == ">=IIB", 1, NA))

TCGA_Test_Dataset$Pathology <- ifelse(TCGA_Test_Dataset$Pathology == "G1", 1,
                                        ifelse(TCGA_Test_Dataset$Pathology == "G2", 2, 3))
TCGA_Test_Dataset$R0.resection <- ifelse(TCGA_Test_Dataset$R0.resection == "Yes", 1,
                                           ifelse(TCGA_Test_Dataset$R0.resection == "No", 0, NA))
TCGA_Test_Dataset$Tumor.site <- ifelse(TCGA_Test_Dataset$Tumor.site == "Body/Tail", 1,
                                         ifelse(TCGA_Test_Dataset$Tumor.site == "Head/Neck", 0, NA))
TCGA_Test_Dataset$AJCC_8th <- ifelse(TCGA_Test_Dataset$AJCC_8th == "<IIB", 0,
                                       ifelse(TCGA_Test_Dataset$AJCC_8th == ">=IIB", 1, NA))

PKUFH_Test_Dataset$Pathology <- ifelse(PKUFH_Test_Dataset$Pathology == "G1", 1,
                                      ifelse(PKUFH_Test_Dataset$Pathology == "G2", 2, 3))
PKUFH_Test_Dataset$R0.resection <- ifelse(PKUFH_Test_Dataset$R0.resection == "Yes", 1,
                                         ifelse(PKUFH_Test_Dataset$R0.resection == "No", 0, NA))
PKUFH_Test_Dataset$Tumor.site <- ifelse(PKUFH_Test_Dataset$Tumor.site == "Body/Tail", 1,
                                       ifelse(PKUFH_Test_Dataset$Tumor.site == "Head/Neck", 0, NA))
PKUFH_Test_Dataset$AJCC_8th <- ifelse(PKUFH_Test_Dataset$AJCC_8th == "<IIB", 0,
                                     ifelse(PKUFH_Test_Dataset$AJCC_8th == ">=IIB", 1, NA))


# Create a list containing training and test sets
list_train_vali_Data <- list(PKUFH_Train_Dataset = PKUFH_Train_Dataset,
                             PKUFH_Test_Dataset = PKUFH_Test_Dataset,
                             TCGA_Test_Dataset = TCGA_Test_Dataset
)
# View list structure
str(list_train_vali_Data)

# Important!!!: Model Construction
library(Mime1)

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$PKUFH_Train_Dataset,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = F,
                       # unicox_p_cutoff = 0.5,
                       candidate_genes = features,
                       mode = 'all',
                       nodesize =5,
                       seed = 123,
                       cores_for_parallel = 10)

pdf("./Result/DFS Plot/DFS_cindex_plot_only_clinical_data.pdf", width = 8, height = 9)  # Adjust width and height as needed

# Run the plotting function
cindex_dis_all(res, 
               validate_set = names(list_train_vali_Data)[-1], 
               order = names(list_train_vali_Data), 
               width = 0.35)

# Close device to save PDF
dev.off()


################### Run DFS model ----------------------
PKUFH_Train_Dataset <- subset(PKUFH_Train_Dataset, select = -c(OS.time,OS))
colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS"] <- "OS"

PKUFH_Test_Dataset <- subset(PKUFH_Test_Dataset, select = -c(OS.time,OS))
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS"] <- "OS"

TCGA_Test_Dataset <- subset(TCGA_Test_Dataset, select = -c(OS.time,OS))
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "DFS.time"] <- "OS.time"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "DFS"] <- "OS"

list_train_vali_Data <- list(PKUFH_Train_Dataset = PKUFH_Train_Dataset,
                             PKUFH_Test_Dataset = PKUFH_Test_Dataset,
                             TCGA_Test_Dataset = TCGA_Test_Dataset
)


library(Mime1)

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$PKUFH_Train_Dataset,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = F,
                       # unicox_p_cutoff = 0.5,
                       candidate_genes = features,
                       mode = 'all',
                       nodesize =5,
                       seed = 123,
                       cores_for_parallel = 10)

pdf("./Result/DFS Plot/DFS_cindex_plot_PCA100_only_clinical_data.pdf", width = 8, height = 15)  # Adjust width and height as needed

# Run the plotting function
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.35)

# Close device to save PDF
dev.off()

##### Pathomics + Clinical Data OS model
PKUFH_Train_Dataset <- subset(PKUFH_Train_Dataset, select = -c(OS.time,OS)) # Clinical data only
PKUFH_Test_Dataset <- subset(PKUFH_Test_Dataset, select = -c(OS.time,OS)) # Clinical data only
TCGA_Test_Dataset <- subset(TCGA_Test_Dataset, select = -c(OS.time,OS)) # Clinical data only

TCGA_Test_pathomics <- read.xlsx("./Processing Result and Dataset/PCA/PCA100/TCGA_data_OS_Model.xlsx")
PKUFH_pathomics <- read.xlsx("./Processing Result and Dataset/PCA/PCA100/PKUFH_data_OS_Model.xlsx")

PKUFH_Train_pathomics <- PKUFH_pathomics[PKUFH_pathomics$ID %in% PKUFH_Train_Dataset$ID,]
# setdiff(TCGA_Test_Dataset$ID, TCGA_Test_pathomics$ID)
PKUFH_Test_pathomics <- PKUFH_pathomics[PKUFH_pathomics$ID %in% PKUFH_Test_Dataset$ID,]

PKUFH_Train_Dataset <- merge(PKUFH_Train_Dataset, PKUFH_Train_pathomics, by.x = "ID", by.y = "ID")
PKUFH_Test_Dataset <- merge(PKUFH_Test_Dataset, PKUFH_Test_pathomics, by.x = "ID", by.y = "ID")
TCGA_Test_Dataset <- merge(TCGA_Test_Dataset, TCGA_Test_pathomics, by.x = "ID", by.y = "ID")
features <- colnames(TCGA_Test_pathomics)
features <- features[! features %in% c("ID", "OS.time", "OS")]
features <- c("ID", "OS.time", "OS", "Pathology","R0.resection","pT.stage","pN.stage","Tumor.site","AJCC_8th",features)
PKUFH_Train_Dataset <- PKUFH_Train_Dataset[,features]
PKUFH_Test_Dataset <- PKUFH_Test_Dataset[,features]
TCGA_Test_Dataset <- TCGA_Test_Dataset[,features]
features <- features[! features %in% c("ID", "OS.time", "OS")]


list_train_vali_Data <- list(PKUFH_Train_Dataset = PKUFH_Train_Dataset,
                             PKUFH_Test_Dataset = PKUFH_Test_Dataset,
                             TCGA_Test_Dataset = TCGA_Test_Dataset
)


library(Mime1)

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$PKUFH_Train_Dataset,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.2,
                       candidate_genes = features,
                       mode = 'all',
                       nodesize =5,
                       seed = 123,
                       cores_for_parallel = 10)

pdf("./Result/OS Plot/OS_cindex_plot_PCA100_with_clinical_data.pdf", width = 8, height = 15)  # Adjust width and height as needed

# Run the plotting function
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.35)

# Close device to save PDF
dev.off()


##### Pathomics + Clinical Data DFS model
PKUFH_Train_Dataset <- subset(PKUFH_Train_Dataset, select = c(ID,Pathology,R0.resection,pT.stage,pN.stage,Tumor.site,AJCC_8th)) # Clinical data only
PKUFH_Test_Dataset <- subset(PKUFH_Test_Dataset, select = c(ID,Pathology,R0.resection,pT.stage,pN.stage,Tumor.site,AJCC_8th)) # Clinical data only
TCGA_Test_Dataset <- subset(TCGA_Test_Dataset, select = c(ID,Pathology,R0.resection,pT.stage,pN.stage,Tumor.site,AJCC_8th)) # Clinical data only

TCGA_Test_pathomics <- read.xlsx("./Processing Result and Dataset/PCA/PCA100/TCGA_data_PFI_Model.xlsx")
PKUFH_pathomics <- read.xlsx("./Processing Result and Dataset/PCA/PCA100/PKUFH_data_DFS_Model.xlsx")

PKUFH_Train_pathomics <- PKUFH_pathomics[PKUFH_pathomics$ID %in% PKUFH_Train_Dataset$ID,]
# setdiff(TCGA_Test_Dataset$ID, TCGA_Test_pathomics$ID)
PKUFH_Test_pathomics <- PKUFH_pathomics[PKUFH_pathomics$ID %in% PKUFH_Test_Dataset$ID,]

PKUFH_Train_Dataset <- merge(PKUFH_Train_Dataset, PKUFH_Train_pathomics, by.x = "ID", by.y = "ID")
PKUFH_Test_Dataset <- merge(PKUFH_Test_Dataset, PKUFH_Test_pathomics, by.x = "ID", by.y = "ID")
TCGA_Test_Dataset <- merge(TCGA_Test_Dataset, TCGA_Test_pathomics, by.x = "ID", by.y = "ID")

colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS"] <- "OS"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS"] <- "OS"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "PFI.time"] <- "OS.time"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "PFI"] <- "OS"


features <- colnames(TCGA_Test_Dataset)
features <- features[! features %in% c("ID", "OS.time", "OS")]
cols <- colnames(TCGA_Test_pathomics)
cols <- cols[!cols %in% c("ID","PFI","PFI.time")]
cols <- c("ID", "OS.time", "OS",features)
PKUFH_Train_Dataset <- PKUFH_Train_Dataset[,cols]
PKUFH_Test_Dataset <- PKUFH_Test_Dataset[,cols]
TCGA_Test_Dataset <- TCGA_Test_Dataset[,cols]
colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "AJCC_8th"] <- "AJCC"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "AJCC_8th"] <- "AJCC"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "AJCC_8th"] <- "AJCC"
features[features == "AJCC_8th"] <- "AJCC"
features



list_train_vali_Data <- list(PKUFH_Train_Dataset = PKUFH_Train_Dataset,
                             PKUFH_Test_Dataset = PKUFH_Test_Dataset,
                             TCGA_Test_Dataset = TCGA_Test_Dataset
)


library(Mime1)

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$PKUFH_Train_Dataset,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.2,
                       candidate_genes = features,
                       mode = 'all',
                       nodesize =5,
                       seed = 123,
                       cores_for_parallel = 10)

pdf("./Result/DFS Plot/DFS_cindex_plot_PCA100_with_clinical_data.pdf", width = 8, height = 15)  # Adjust width and height as needed

# Run the plotting function
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.35)

# Close device to save PDF
dev.off()
