##############################
##### required libraries #####
##############################
# library(data.table)
# library(randomForest)
# library(caret)
# library(PRROC)
# library(caret)
# library(ada)
# library(ROCR)
# library(mlbench)
# library(caret)
# library(tidyverse)

#############################
#### required functions #######
#############################
# create_train_test <- function(data, size = 0.8, train = TRUE) {
#   n_row = nrow(data)
#   total_row = size * n_row
#   train_sample <- 1: total_row
#   if (train == TRUE) {
#     return (data[train_sample, ])
#   } else {
#     return (data[-train_sample, ])
#   }
# }

########################################
######## TRAINING THE MODEL ############
########################################

data <- read.table("IPO323.reference.metrics.table", sep="\t", h=T)
data <- data[,-1:-3]
# ignore SNPs in regions of high TE coverage
data$snp <- with(data,ifelse(te_cov>0.8,"NA",snp))
data$snp <- as.numeric(as.character(data$snp))

# define translocations dataset
tra_data <- subset(data, colnames %in% c("DHH", "DTA", "DTB", "DTC" ,"DTH" ,"DTM" ,"DTP", "DTT", "DTX" ,"DXX" ,"DYC", "RII", "RIJ" ,"RIX" ,"RLB" ,"RLC" ,"RLG", "RLX" ,"RSX" ,"RXX" ,"RYN", "XXX", "nested","gene","TE","GC", "H3K9", "H3K4", "H3K27", "SNP", "cM","core" ,"accessory","singleton", "effector", "TRA"))
# define translocations as factors for the model
tra_data$TRA <- as.factor(tra_data$TRA)
# order the dataset for consistency
tra_data <- tra_data[,c("DHH", "DTA", "DTB", "DTC" ,"DTH" ,"DTM" ,"DTP", "DTT", "DTX" ,"DXX" ,"DYC", "RII", "RIJ" ,"RIX" ,"RLB" ,"RLC" ,"RLG", "RLX" ,"RSX" ,"RXX" ,"RYN", "XXX", "nested","gene","TE","GC","cM","core" ,"accessory", "effector", "TRA")]
# define indels dataset
indel_data <- subset(data, colnames %in% c("DHH", "DTA", "DTB", "DTC" ,"DTH" ,"DTM" ,"DTP", "DTT", "DTX" ,"DXX" ,"DYC", "RII", "RIJ" ,"RIX" ,"RLB" ,"RLC" ,"RLG", "RLX" ,"RSX" ,"RXX" ,"RYN", "XXX", "nested","gene","TE","GC", "H3K9", "H3K4", "H3K27", "SNP", "cM","core" ,"accessory","singleton", "effector", "INDEL"))
# define indels as factors for the model
indel_data$INDEL <- as.factor(indel_data$INDEL)
# order the dataset for consistency
indel_data <- indel_data[,c("DHH", "DTA", "DTB", "DTC" ,"DTH" ,"DTM" ,"DTP", "DTT", "DTX" ,"DXX" ,"DYC", "RII", "RIJ" ,"RIX" ,"RLB" ,"RLC" ,"RLG", "RLX" ,"RSX" ,"RXX" ,"RYN", "XXX", "nested","gene","TE","GC", "cM","core" ,"accessory", "effector", "INDEL")]
# set seed for reproducibility
# set.seed(123)
set.seed(1234)
# define train and test datasets for both indels and translocations
tra_train <- create_train_test(tra_data, 0.8, train = TRUE)
tra_test <- create_train_test(tra_data, 0.2, train = FALSE)
indel_train <- create_train_test(indel_data, 0.8, train = TRUE)
indel_test <- create_train_test(indel_data, 0.2, train = FALSE)
# define training control parameter
train.control <- trainControl(method = "repeatedcv", number = 10, repeats=3)

### RF METHOD
# set tuning parameters
tune1 <- data.frame("mtry"=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
# train the random forest model for indels and translocations indepedently
model1 <- train(TRA ~ ., data = tra_train, method = "rf", metric="Accuracy", num.threads = 4, tuneGrid = tune1, trControl = train.control)
#saveRDS(model1, "RF_model_tra.rds")
model2 <- train(INDEL ~ ., data = indel_train, method = "rf", metric="Accuracy", num.threads = 4, tuneGrid = tune1, trControl = train.control)
#saveRDS(model2, "RF_model_indel.rds")

### GLM METHOD
model3 <- train(TRA ~ ., data = tra_train, method = "glm", metric="Accuracy", family = "binomial", trControl = train.control)
#saveRDS(model3, "GLM_model_tra.rds")
model4 <- train(INDEL ~ ., data = indel_train, method = "glm", metric="Accuracy", family = "binomial", trControl = train.control)
#saveRDS(model4, "GLM_model_indel.rds")

### ADA METHOD
# set tuning parameters
tune2 <- data.frame("iter"=c(100, 100, 100, 1000, 1000, 1000, 3000, 3000, 3000), "maxdepth"=c(1,5,20,1,5,20,1,5,20), "nu"=c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01))
model5 <- train(TRA ~ ., data = tra_train, method = "ada", metric="Accuracy", num.threads = 4, tuneGrid = tune2, trControl = train.control)
#saveRDS(model5, "ADA_model_tra.rds")
model6 <- train(INDEL ~ ., data = indel_train, method = "ada", metric="Accuracy", num.threads = 4, tuneGrid = tune2, trControl = train.control)
#saveRDS(model6, "ADA_model_indel.rds")

### GBM METHOD
# set tuning parameters
tune3 <- data.frame("interaction.depth"=c(1, 5, 10, 1, 5, 10, 1, 5, 10, 1, 5, 10, 1, 5, 10, 1, 5, 10), "n.trees"=c(50, 50, 50, 500, 500, 500, 1000, 1000, 1000, 50, 50, 50, 500, 500, 500, 1000, 1000, 1000), "shrinkage"=c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01), "n.minobsinnode"=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5))
model7 <- train(TRA ~ ., data = tra_train, method = "gbm", metric="Accuracy", num.threads = 4, tuneGrid = tune3, trControl = train.control)
#saveRDS(model7, "GBM_model_tra.rds")
model8 <- train(INDEL ~ ., data = indel_train, method = "gbm", metric="Accuracy", num.threads = 4, tuneGrid = tune3, trControl = train.control)
#saveRDS(model8, "GBM_model_indel.rds")

### estimate variables importance for each model
varImp(model1, scale = F)$importance

######################################################
### evaluate model performance on the test dataset ###
######################################################
tra_test <- predict(model1, tra_test, type = "prob")
indel_test <- predict(model2, indel_test, type = "prob")
# use ROCR's prediction function to format the data
tra_pred <- prediction(tra_test[,2], tra_test$TRA)
indel_pred <- prediction(indel_test[,2], indel_test$INDEL)
# uses the prediction object to derive performance metrics (in this case tpr and fpr)
perf1 <- performance(tra_pred, "tpr", "fpr")
fg1 <- tra_test[, 2]
bg1 <- tra_test[, 1]
perf2 <- performance(indel_pred, "tpr", "fpr")
fg2 <- indel_test[, 2]
bg2 <- indel_test[, 1]
# compute the PR Curve
pr1 <- pr.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)
pr2 <- pr.curve(scores.class0 = fg2, scores.class1 = bg2, curve = T)
# plot the ROC and PROC curves
plot(perf1, col="orange", lty = 2, lwd = 2)
plot(pr1, col="orange", lty = 2, lwd = 2, main="", auc.main=F)
# predicting indels and translocations given the test data
tra_predValid <- predict(model1, tra_test, type = "raw")
indel_predValid <- predict(model2, indel_test, type = "raw")
# compute the confusion matrix
tra_result <- confusionMatrix(tra_predValid, tra_test$TRA, positive = "1")
indel_result <- confusionMatrix(indel_predValid, indel_test$INDEL, positive = "1")
# extract each model summary statistics
stats1 <- data.frame(cbind("TRA_RF",t(tra_result$overall), t(tra_result$byClass)))
table1 <- as.table(tra_result1)
stats2 <- data.frame(cbind("INDEL_RF",t(indel_result$overall), t(indel_result$byClass)))
table2 <- as.table(indel_result2)

