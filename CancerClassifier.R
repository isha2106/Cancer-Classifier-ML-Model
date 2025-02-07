# Load the data
data <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vSx3gPbBa0-_s2mPqlXij2sTVdoF5hcp7zb5NIAnaeMuwlLZ5PXmijMaOWSk0nVDw/pub?gid=1571246026&single=true&output=csv")

# Load the labels
label_data <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQvYplBm1JbPw-VSFB79N01SeR4lWTCh9yl3WUiFtiqHVenbYTijM5ryqJ5ELzAew/pub?gid=1180405259&single=true&output=csv")

# Merge the data with the labels
cancer_data <- merge(data, label_data, by = "X")

if (!require("dplyr", character.only = TRUE)) {
  install.packages("dplyr", dependencies = TRUE)
}
library(dplyr)

# Dropping column X (samplename)
cancer_data <- select(cancer_data, -X)

if (!require("ggplot2", character.only = TRUE)) {
  install.packages("ggplot2", dependencies = TRUE)
}
library(ggplot2)

# Plot a histogram
ggplot(cancer_data, aes(x = Class)) + 
  geom_bar(fill = "steelblue") + 
  theme_minimal() + 
  labs(title = "Number of cases for different types of Cancer", x = "Class", y = "Count") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Checking for missing values
sum(is.na(cancer_data))

# Impute missing values with mean of each column
cancer_data[] <- lapply(cancer_data, function(x) {
  if(is.numeric(x)) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
  }
  return(x)
})

# Check again for missing values
sum(is.na(cancer_data))

# Function for finding outliers
outliers <- function(data, column) {
  
  # Calculate mean and standard deviation
  mu <- mean(data[[column]], na.rm = TRUE)
  sd <- sd(data[[column]], na.rm = TRUE)
  
  # Calculate z-score
  z_scores <- (data[[column]] - mu) / sd
  
  # Identify z-scores above the threshold
  outlier_threshold <- which(abs(z_scores) > 2.5)
  
  return(outlier_threshold)
}

# Call the function for each gene
cols <- names(cancer_data)[-c(1, ncol(cancer_data))]
outliers_list <- list()

for(column in cols) {
  outliers_list[[column]] <- outliers(cancer_data, column)
}

total_outliers <- sum(lengths(outliers_list))
print(paste("The total number of outliers are:", total_outliers))

# Visualize the outlier counts
outlier_counts <- sapply(outliers_list, length)
barplot(outlier_counts, main = "Outlier Counts per Gene", xlab = "Genes", ylab = "Number of Outliers", las = 2)

# Let's visualize genes with outlier counts above a certain threshold, e.g., greater than 10
high_outlier_genes <- names(which(outlier_counts > 20))
if (length(high_outlier_genes) > 0) {
  boxplot(cancer_data[, high_outlier_genes], las = 2, main = "Boxplots for Genes with High Outlier Counts", xlab = "Genes", ylab = "Expression Levels")
} else {
  print("No genes have more than 10 outliers.")
}

if (!require("reshape2", character.only = TRUE)) {
  install.packages("reshape2", dependencies = TRUE)
}
library(reshape2)

# One hot encoding
cancer_data$Class <- as.factor(cancer_data$Class)

# Create dummy variables
dummy_class <- model.matrix(~ Class - 1, data = cancer_data)

# Convert all data to numeric
gene_data <- select(cancer_data, -Class)
gene_data <- as.data.frame(lapply(gene_data, as.numeric))

# Calculate correlations
cor_analysis <- (cor(cbind(gene_data, dummy_class), use = "complete.obs"))
cor_genes <- ncol(gene_data)
cor_class <- ncol(dummy_class)
final_cor <- cor_analysis[1:cor_genes, (cor_genes+1):(cor_genes+cor_class)]

# Melt and format the matrix 
melt_cor <- melt(final_cor, varnames = c("Gene", "Class"))
melt_cor$Class <- gsub("Class", "", melt_cor$Class)

# Find the least 10 correlated genes for each class
order_data <- melt_cor[order(melt_cor$Class, melt_cor$value), ]
group_data <- group_by(order_data, Class)
least_cor_genes <- summarise(group_data, least_genes = Gene[1:5], least_val = value[1:5], .groups = 'drop')

print(least_cor_genes)

# Find common least correlated genes between the classes
genes_per_class <- by(least_cor_genes, least_cor_genes$Class, function(df) {
  as.character(df$least_genes)
})
common_genes <- Reduce(intersect, genes_per_class)
print(common_genes)


if (!require("umap", character.only = TRUE)) {
  install.packages("umap", dependencies = TRUE)
}
library(umap)

set.seed(3434)

# UMAP dimensionality reduction
umap_plot <- umap(gene_data)

# Convert for ggplot
umap_df <- as.data.frame(umap_plot$layout)
umap_df$Class <- cancer_data$Class

# Plot a UMAP
ggplot(umap_df, aes(x = V1, y = V2, color = Class)) +
  geom_point(alpha = 0.7) +
  labs(title = "UMAP for Correlation between Genes of each Class", x = "Genes", y = "Genes") +
  theme_minimal() +
  theme(legend.position = "right")


# Check correlation between genes
gene_cor_matrix <- cor(gene_data, use = "complete.obs")
gene_cor_matrix <- melt(gene_cor_matrix)

ggplot(gene_cor_matrix, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Correlation Values", x = "Correlation Value", y = "Frequency") +
  theme_minimal()

# Perform log transform
log_tranform_data <- select(cancer_data, -Class)
log_tranform_data <- log2(log_tranform_data + 1)
log_tranform_data$Class <- cancer_data$Class

summary(log_tranform_data)

# Function for normalization
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) }

# Call the function
norm_data <- select(log_tranform_data, -Class)
norm_data <- as.data.frame(lapply(norm_data, normalize))
norm_data$Class <- log_tranform_data$Class

summary(norm_data)

# Function for finding outliers
outliers <- function(data, column) {
  
  # Calculate mean and standard deviation
  mu <- mean(data[[column]], na.rm = TRUE)
  sd <- sd(data[[column]], na.rm = TRUE)
  
  # Calculate z-score
  z_scores <- (data[[column]] - mu) / sd
  
  # Identify z-scores above the threshold
  outlier_threshold <- which(abs(z_scores) > 2.5)
  
  return(outlier_threshold)
}

# Call the function for each gene
cols <- names(norm_data)[-c(1, ncol(norm_data))]
outliers_list <- list()

for(column in cols) {
  outliers_list[[column]] <- outliers(norm_data, column)
}
#outliers_list

# Total number of outliers
total_outliers <- sum(lengths(outliers_list))
print(paste("The total number of outliers are:", total_outliers))

# Remove columns with only 0 values
normalized_data <- select(norm_data, -Class)
zero_cols <- sapply(normalized_data, function(x) length(unique(x)) == 1)
constant_column_names <- names(normalized_data)[zero_cols]
normalized_data <- normalized_data[, !zero_cols]

print(paste(constant_column_names, "removed because they only have zero values"))

# Apply PCA
pca_results <- prcomp(normalized_data)
str(pca_results)
summary(pca_results)

# Plot PCA graph
plot(pca_results, type = "l")

# Making a new dataset
pc_scores <- as.data.frame(pca_results$x)
pc_scores$Class <- norm_data$Class

# Components that reach 90% variance
pc_comp <- cumsum(pca_results$sdev^2 / sum(pca_results$sdev^2))
num_components <- which(pc_comp >= 0.90)[1]   
pca_dataset <- pc_scores[, 1:num_components]
pca_dataset$Class <- norm_data$Class 

# Normalized Dataset
normalized_data$Class <- norm_data$Class
str(normalized_data)

# PCA Dataset
str(pca_dataset)

## Model Construction
if (!require("caret", character.only = TRUE)) {
  install.packages("caret", dependencies = TRUE)
}
library(caret)

set.seed(3434)
# Subsets based on index and classes
norm_train_index <- createDataPartition(normalized_data$Class, p = 0.8, list = FALSE)

# Training set
norm_train_subset <- normalized_data[norm_train_index, ]

# Testing set
norm_test_subset <- normalized_data[-norm_train_index, ]

set.seed(3434)
# Subsets based on index and classes
pca_train_index <- createDataPartition(pca_dataset$Class, p = 0.8, list = FALSE)

# Training set
pca_train_subset <- pca_dataset[pca_train_index, ]

# Testing set
pca_test_subset <- pca_dataset[-pca_train_index, ]

if (!require("kernlab", character.only = TRUE)) {
  install.packages("kernlab", dependencies = TRUE)
}
library(kernlab)

# SVM model for norm_train_subset
set.seed(3434)
svm_model_norm <- ksvm(Class ~ ., data = norm_train_subset, kernel = "rbfdot")

# SVM model for pca_train_subset
set.seed(3434)
svm_model_pca <- ksvm(Class ~ ., data = pca_train_subset, kernel = "rbfdot")

if (!require("gbm", character.only = TRUE)) {
  install.packages("gbm", dependencies = TRUE)
}
library(gbm)

# Gradient Boosting model for normalized_data
set.seed(3434)
gbm_model_norm <- train(Class ~ ., data = norm_train_subset, method = "gbm", verbose = FALSE)

# Gradient Boosting model for normalized_data
set.seed(3434)
gbm_model_pca <- train(Class ~ ., data = pca_train_subset, method = "gbm", verbose = FALSE)

if (!require("randomForest", character.only = TRUE)) {
  install.packages("randomForest", dependencies = TRUE)
}
library(randomForest)

# Random Forest model for normalized_data
set.seed(3434)
rf_model_norm <- randomForest(Class ~ ., data = norm_train_subset, verbose = FALSE)

# Random Forest model for pca_train_subset
set.seed(3434)
rf_model_pca <- randomForest(Class ~ ., data = pca_train_subset, verbose = FALSE)

importance(rf_model_norm)

# Plot variable importance
varImpPlot(rf_model_norm)
varImpPlot(rf_model_pca)

# Evaluate SVM model for norm_test_subset
set.seed(3434)
svm_predict_norm_test <- predict(svm_model_norm, norm_test_subset)
confusionMatrix(svm_predict_norm_test, norm_test_subset$Class)

# Evaluate SVM model for pca_test_subset
set.seed(3434)
svm_predict_pca_test <- predict(svm_model_pca, pca_test_subset)
confusionMatrix(svm_predict_pca_test, pca_test_subset$Class)

# Evaluate Gradient Boosting model for norm_test_subset
set.seed(3434)
gb_predict_norm_test <- predict(gbm_model_norm, norm_test_subset)
confusionMatrix(gb_predict_norm_test, norm_test_subset$Class)

# Evaluate Gradient Boosting model for pca_test_subset
set.seed(3434)
gb_predict_pca_test <- predict(gbm_model_pca, pca_test_subset)
confusionMatrix(gb_predict_pca_test, pca_test_subset$Class)

# Evaluate Random Forest model for norm_test_subset
set.seed(3434)
rf_predict_norm_test <- predict(rf_model_norm, norm_test_subset)
confusionMatrix(rf_predict_norm_test, norm_test_subset$Class)

# Evaluate Random Forest model for pca_test_subset
set.seed(3434)
rf_predict_pca_test <- predict(rf_model_pca, pca_test_subset)
confusionMatrix(rf_predict_pca_test, pca_test_subset$Class)

# Creating train control
set.seed(3434)
train_control <- trainControl(method = "cv", number = 5, verboseIter = FALSE, savePredictions = "final",                               classProbs = TRUE, summaryFunction = multiClassSummary)

if (!require("caret", character.only = TRUE)) {
  install.packages("caret", dependencies = TRUE)
}
library(caret)

# SVM grid
svm_grid <- expand.grid(C = c(0.1, 1, 10), sigma = c(0.001, 0.01, 0.05))

# k-fold cross-validation SVM model for normalized_data
set.seed(3434)
cross_svm_model_norm <- train(Class ~ ., data = normalized_data, method = "svmRadial", trControl =       
                                train_control, metric = "Accuracy", tuneGrid = svm_grid)

# k-fold cross-validation SVM model for pca_dataset
set.seed(3434)
cross_svm_model_pca <- train(Class ~ ., data = pca_dataset, method = "svmRadial", trControl = 
                               train_control, metric = "Accuracy", tuneGrid = svm_grid)

# Gradient Boosting grid
gbm_grid <- expand.grid(n.trees = c(50, 100), interaction.depth = c(1, 3), shrinkage = c(0.1),  n.minobsinnode = c(5, 10))

# k-fold cross-validation Gradient Boosting model for normalized_data
set.seed(3434)
cross_gbm_model_norm <- train(Class ~ ., data = normalized_data, method = "gbm", trControl = 
                                train_control, metric = "Accuracy", verbose = FALSE, tuneGrid = gbm_grid)

# k-fold cross-validation Gradient Boosting model for pca_dataset
set.seed(3434)
cross_gbm_model_pca <- train(Class ~ ., data = pca_dataset, method = "gbm", trControl = train_control,   
                             metric = "Accuracy", verbose = FALSE, tuneGrid = gbm_grid)

# Random forest grid
rf_grid <- expand.grid(mtry = c(2, 4, 6))

# k-fold cross-validation Gradient Boosting model for normalized_data
set.seed(3434)
cross_rf_model_norm <- train(Class ~ ., data = normalized_data, method = "rf", trControl = train_control, 
                             metric = "Accuracy", verbose = FALSE, tuneGrid = rf_grid)

# k-fold cross-validation Gradient Boosting model for pca_dataset
set.seed(3434)
cross_rf_model_pca <- train(Class ~ ., data = pca_dataset, method = "rf", trControl = train_control,   
                            metric = "Accuracy", verbose = FALSE, tuneGrid = rf_grid)

create_ensemble <- function(rf_model, gbm_model, svm_model, test_data) {
  
  # Predictions for Random Forest
  rf_predictions <- predict(rf_model, test_data)
  
  # Predictions for Gradient Boosting
  gbm_predictions <- predict(gbm_model, newdata = test_data)
  
  # Predictions for SVM
  svm_predictions <- predict(svm_model, test_data)
  
  # Combine predictions into a matrix
  combined_predictions <- cbind(as.character(rf_predictions),
                                as.character(gbm_predictions),
                                as.character(svm_predictions))
  
  # Perform majority vote with random selection in case of a tie
  ensemble_predictions <- apply(combined_predictions, 1, function(row) {
    freq <- table(row)
    return(names(which.max(freq)))
  })
  return(ensemble_predictions)
}

norm_ensemble <- create_ensemble(cross_rf_model_norm, cross_gbm_model_norm, cross_svm_model_norm, norm_test_subset)

pca_ensemble <- create_ensemble(cross_rf_model_pca, cross_gbm_model_pca, cross_svm_model_pca, pca_test_subset)

norm_ensemble_accuracy <- sum(norm_ensemble == norm_test_subset$Class) / length(norm_test_subset$Class) 
print(paste("Ensemble Accuracy:", norm_ensemble_accuracy))

pca_ensemble_accuracy <- sum(pca_ensemble == pca_test_subset$Class) / length(pca_test_subset$Class) 
print(paste("Ensemble Accuracy:", pca_ensemble_accuracy))



























