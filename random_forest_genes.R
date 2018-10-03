## ------Error rate for genes-------#
#setwd("E:\\RShinyapp\\networksanalysis\\code")
#setwd("C://Networks//code//")
require("randomForest")
require(igraph)
require(BNN)
print("check 1")
diff_exp_genes <- read.csv("Diff_Exp_Data.csv")
diff_exp_genes <- diff_exp_genes[,-1]
feat_select_genes <- read.csv("Feature_Selection.csv")
feat_select_genes <- feat_select_genes[,-1]
drop_cols <- c("IDENTIFIER")
feat_select_genes <- feat_select_genes[,!(names(feat_select_genes) %in% drop_cols)]
expression_data_gds5037 <- read.csv("exprsdata.csv")
drops <- c("X")
expression_df_gds5037 <- expression_data_gds5037[ , !(names(expression_data_gds5037) %in% drops)]
print("check 2")
metadata_gds5037 <- read.csv("metadata.csv")
metadata_gds5037 <- metadata_gds5037[,-1]
diseased <- as.data.frame(read.csv("Diseased.csv"))
rf_diff_genes <- list()
rf_final_genes <- list()
rf_mod_j <- list()
rf_result <- list()
print("check 3")
for(i in 1:nrow(diff_exp_genes)){
  rf_diff_genes[[i]] <- diff_exp_genes[i,ncol(diff_exp_genes)]
  rf_diff_genes[[i]] <- as.data.frame(rf_diff_genes[[i]])
  colnames(rf_diff_genes[[i]]) <- 'ID_REF'
  get_exp_data <- merge(expression_df_gds5037,rf_diff_genes[[i]], by = c('ID_REF'))
  drop_cols <- c("ID_REF", "IDENTIFIER")
  exp_df_for_rf <- get_exp_data[,!(names(get_exp_data) %in% drop_cols)]
  exp_df_for_rf_t <- t(exp_df_for_rf)
  exp_df_class <- cbind(metadata_gds5037$disease.state, exp_df_for_rf_t)
  exp_df_class <- as.data.frame(exp_df_class)
  exp_df_class <- t(exp_df_class)
  d <- as.matrix(diseased$V1)
  if(e$algo == 'bor'){
    print(exp_df_class[1,])
    print(d)
    rf <- randomForest(d~ exp_df_class[1,], ntree=1000,importance=T)
    rf_mod_j[[i]] <- rf$mse[1000]
  }
  if(e$algo == 'bnn'){
    bnn <- BNNsel(d,exp_df_class[1,],total_iteration = 10000)
    b <- 0
    for(j in (1:nrow(d))){
      if(j<=length(bnn$fit)){
        b <- b+ (1-bnn$fit[j])
      }
    }
    rf_mod_j[[i]] <- b/length(bnn$fit)
  }
}
print("check 4")
rf_result[[1]] <- rf_mod_j
rf_mod_j <- list()
for(i in 1:nrow(feat_select_genes)){
  rf_final_genes[[i]] <- feat_select_genes[i,]
  rf_final_genes[[i]] <- as.data.frame(rf_final_genes[[i]])
  colnames(rf_final_genes[[i]]) <- 'ID_REF'
  get_exp_data <- rf_final_genes[[i]]
  drop_cols <- c("ID_REF", "IDENTIFIER")
  exp_df_for_rf <- get_exp_data[,!(names(get_exp_data) %in% drop_cols)]
  exp_df_for_rf_t <- t(exp_df_for_rf)
  exp_df_class <- cbind(metadata_gds5037$disease.state, exp_df_for_rf_t)
  exp_df_class <- as.data.frame(exp_df_class)
  d <- as.matrix(diseased$V1)
  if(e$algo == 'bor'){
    rf <- randomForest(d~ exp_df_class[,1], ntree=1000,importance=T)
    rf_mod_j[[i]] <- rf$mse[1000]
  }
  if(e$algo == 'bnn'){
    bnn <- BNNsel(d,exp_df_class[,1],total_iteration = 10000)
    b <- 0
    for(j in (1:nrow(d))){
      if(j<=length(bnn$fit)){
        b <- b+ (1-bnn$fit[j])
      }
    }
    rf_mod_j[[i]] <- b/length(bnn$fit)
  }
}

rf_result[[2]] <- rf_mod_j
print(rf_result[[2]])
if(is.null(unlist(rf_result[[1]]))){
  m$V1 <- data.frame(matrix(0, ncol = 1, nrow = nrow(rf_result[[2]])))
}else{
  m$V1 <- data.frame(matrix(unlist(rf_result[[1]]), byrow=T))
}
if(is.null(unlist(rf_result[[2]]))){
  m$V2 <- data.frame(matrix(0, ncol = 1, nrow = nrow(rf_result[[1]])))
}else{
  m$V2 <- data.frame(matrix(unlist(rf_result[[2]]), byrow=T))
}

rf_module_df_temp <- as.data.frame(do.call(cbind, m))
rf_module_dfs <- apply(rf_module_df_temp,2,as.character)
write.csv(rf_module_dfs[,-1], "rf_genes.csv", row.names = F)
rf_csv <- read.csv("rf_genes.csv", header = T)
