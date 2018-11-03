## ------Random FOrest error rate for each modules for comaprison-------#
#setwd("E:\\RShinyapp\\networksanalysis\\code")
#setwd("C://Networks//code//")
require(igraph)
edgelist <- read.csv("Edgelist_top_10prcnt_genes.csv")
# in case of correlation converting to absolute values
edgelist["abs_value"] <- abs(edgelist[,3])
edgelist = edgelist[,-3]


graph_obj <- graph.data.frame(edgelist)
imap_comm_det <- cluster_infomap(graph_obj, e.weights = edgelist[,3], v.weights = edgelist[,1], nb.trials = 10,
                                 modularity = F)
#total number of communities
total_com <- max(imap_comm_det$membership)
#predictive analysis per module
expression_data_gds5037 <- read.csv("exprsdata.csv")

#NOTE :  WHEN taking edata, make the following changes -
#expression_data_gds5037 <- edata

# All this X, ID for gds5037
drops <- c("X")
expression_df_gds5037 <- expression_data_gds5037[ , !(names(expression_data_gds5037) %in% drops)]
#metadata_gds5037 <- mdata
# when reading directly use col -1
metadata_gds5037 <- read.csv("metadata.csv")
metadata_gds5037 <- metadata_gds5037[,-1]
diseased <- as.data.frame(read.csv("Diseased.csv"))
require("randomForest")

rf_mod_i <- list()
rf_mod_j <- list()
rf_module_df <- list()
for(i in (1:total_com))
{
  module_i <- as.data.frame(imap_comm_det[i])
  if(nrow(module_i) > 0){
    print(i)
    colnames(module_i) <- 'ID_REF'
    for(j in 1:nrow(module_i)){
      rf_module_df[[j]] <- module_i[j,]
      rf_module_df[[j]] <- as.data.frame(rf_module_df[[j]])
      colnames(rf_module_df[[j]]) <- 'ID_REF'
      get_exp_data <- merge(expression_df_gds5037,rf_module_df[[j]], by = c('ID_REF'))
      drop_cols <- c("ID_REF", "IDENTIFIER")
      exp_df_for_rf <- get_exp_data[,!(names(get_exp_data) %in% drop_cols)]
      exp_df_for_rf_t <- t(exp_df_for_rf)
      exp_df_class <- cbind(metadata_gds5037$disease.state, exp_df_for_rf_t)
      exp_df_class <- as.data.frame(exp_df_class)
      d <- as.matrix(diseased$V1)
      rf <- randomForest(d~. , data=exp_df_class, ntree=1000,importance=T)
      rf_mod_j[[j]] <- rf$mse[1000]

    }
    rf_mod_i[[i]] <- rf_mod_j
  }
}

rf_module_df_temp <- as.data.frame(do.call(cbind, rf_mod_i))
rf_module_dfs <- apply(rf_module_df_temp,2,as.character)
write.csv(rf_module_dfs, "rf_modules.csv", row.names = F)
rf_csv <- read.csv("rf_modules.csv", header = T)

boxplot(rf_csv, xlab = "Modules", ylab = "Random forest error rate")
