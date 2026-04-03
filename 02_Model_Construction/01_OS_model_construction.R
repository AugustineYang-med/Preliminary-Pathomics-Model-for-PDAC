set.seed(123)  # Set random seed for reproducibility
library(caret)

OS_normalized_PKUFH_data$OS <- as.factor(OS_normalized_PKUFH_data$OS) # Ensure OS is a factor (for 0/1 events)
class(OS_normalized_PKUFH_data$OS)  # Should be "factor"
table(OS_normalized_PKUFH_data$OS)  # Check event distribution

identical(colnames(OS_normalized_PKUFH_data), colnames(OS_normalized_TCGA_data))

# Generate all combinations from AAA to ZZZ (26^3 = 17576 possibilities)
all_combinations <- apply(
  expand.grid(LETTERS, LETTERS, LETTERS), 
  1, 
  paste, 
  collapse = ""
)
three_letter_names <- all_combinations[1:3330] # Take first 3330 (in alphabetical order)

# Get total number of columns
n_cols <- ncol(OS_normalized_PKUFH_data)
# Replace column names of last 3330 columns
start_col <- n_cols - 3330 + 1  # Calculate starting column position
colnames(OS_normalized_PKUFH_data)[start_col:n_cols] <- three_letter_names
colnames(OS_normalized_TCGA_data)[start_col:n_cols] <- three_letter_names
identical(colnames(OS_normalized_PKUFH_data), colnames(OS_normalized_TCGA_data))

pathomics_features <- colnames(OS_normalized_TCGA_data)
pathomics_features <- pathomics_features[pathomics_features != c("ID", "OS.time", "OS")]


# PKUFH dataset stratified sampling (based on OS event)
PKUFH_train_index <- createDataPartition(y = OS_normalized_PKUFH_data$OS,  # Must be a factor or numeric vector
                                   p = 0.7,                          # Training set ratio (70%)
                                   list = TRUE,                      # Return list format
                                   times = 1                         # Number of samplings (1)
                                   )

# Extract training and test sets
PKUFH_train_data <- OS_normalized_PKUFH_data[PKUFH_train_index$Resample1, ]  # Training set (70%)
PKUFH_test_data  <- OS_normalized_PKUFH_data[-PKUFH_train_index$Resample1, ] # Test set (30%)
identical(colnames(PKUFH_train_data), colnames(PKUFH_test_data))


# Create a list containing training and test sets
list_train_vali_Data <- list(PKUFH_Train_Dataset = PKUFH_train_data,
                             PKUFH_Test_Dataset = PKUFH_test_data,
                             TCGA_Validation_Dataset = OS_normalized_TCGA_data
                             )

# View list structure
str(list_train_vali_Data)

# Check event proportions in training and test sets
prop.table(table(list_train_vali_Data$PKUFH_Train_Dataset$OS))  # Training set proportion
prop.table(table(list_train_vali_Data$PKUFH_Test_Dataset$OS))   # Test set proportion


# Important!!!: Model Construction
library(Mime1)

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$PKUFH_Train_Dataset,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.001,
                       candidate_genes = pathomics_features,
                       mode = 'all',
                       nodesize =5,
                       seed = 123,
                       cores_for_parallel = 10)

pdf("cindex_plot.pdf", width = 8, height = 15)  # Adjust width and height as needed

# Run the plotting function
cindex_dis_all(res, 
               validate_set = names(list_train_vali_Data)[-1], 
               order = names(list_train_vali_Data), 
               width = 0.35)

# Close device to save PDF
dev.off()

#####################################################------------------------------
######################### Loop function to generate models with different parameters
#####################################################------------------------------


# Define parameter ranges
nodesize_range <- 5:10
unicox_p_cutoff_values <- c(0.01, 0.001)

# Create empty list to store results
results <- list()

# Loop execution
for (nodesize in nodesize_range) {
  for (unicox_p_cutoff in unicox_p_cutoff_values) {
    # Print current parameter combination
    cat("Running with nodesize =", nodesize, "and unicox_p_cutoff =", unicox_p_cutoff, "\n")
    
    # Execute function
    current_result <- ML.Dev.Prog.Sig(
      train_data = list_train_vali_Data$PKUFH_Train_Dataset,
      list_train_vali_Data = list_train_vali_Data,
      unicox.filter.for.candi = TRUE,
      unicox_p_cutoff = unicox_p_cutoff,
      candidate_genes = pathomics_features,
      mode = 'all',
      nodesize = nodesize,
      seed = 123
    )
    
    # Store results using parameter combination as name
    result_name <- paste0("nodesize_", nodesize, "_cutoff_", unicox_p_cutoff)
    results[[result_name]] <- current_result
  }
}

# Return all results
results

