set.seed(123)  # Set random seed for reproducibility
library(openxlsx)
PKUFH_Train_Dataset <- read.xlsx("./Result/PCA/PCA100/PKUFH_train_data_DFS_Model.xlsx")
PKUFH_Test_Dataset <- read.xlsx("./Result/PCA/PCA100/PKUFH_test_data_DFS_Model.xlsx")
TCGA_Test_Dataset <- read.xlsx("./Data/TCGA-PAAD_pathomics_survival_data.xlsx")
colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Train_Dataset)[colnames(PKUFH_Train_Dataset) == "DFS"] <- "OS"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS.time"] <- "OS.time"
colnames(PKUFH_Test_Dataset)[colnames(PKUFH_Test_Dataset) == "DFS"] <- "OS"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "PFI.time"] <- "OS.time"
colnames(TCGA_Test_Dataset)[colnames(TCGA_Test_Dataset) == "PFI"] <- "OS"

cols <- colnames(TCGA_Test_Dataset)
cols <- cols[cols != c("ID", "OS.time", "OS")]

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
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = cols,
                       mode = 'all',
                       nodesize =5,
                       seed = 123,
                       cores_for_parallel = 10)

pdf("./Result/DFS Plot/DFS_cindex_plot_PCA100.pdf", width = 8, height = 15)  # Adjust width and height as needed

# Run the plotting function
cindex_dis_all(res, 
               validate_set = names(list_train_vali_Data)[-1], 
               order = names(list_train_vali_Data), 
               width = 0.35)

# Close device to save PDF
dev.off()