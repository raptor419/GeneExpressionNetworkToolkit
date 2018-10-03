## ------INFO map community detection-------
#setwd("E:\\RShinyapp\\networksanalysis\\code")
#setwd("C://Networks//code//")
require(igraph)
edgelist <- read.csv("Edgelist_top_10prcnt_genes.csv")
graph_obj <- graph.data.frame(edgelist)
imap_comm_det <- cluster_infomap(graph_obj, e.weights = edgelist[,3], v.weights = edgelist[,1], nb.trials = 10,
                modularity = TRUE)
#total number of communities
total_com <- max(imap_comm_det$membership)
#class(total_com)
#sprintf("total number of modules = %d", total_com)
#cat("\n")
#cat("total \n", total_comm)
#sprintf("Modularity = %f", imap_comm_det$modularity)
#cat("\n")
#view each community
#imap_comm_det[1]
#imap_comm_det
#print(imap_comm_det[as.numeric(input$module_num)])

#edgelist[(as.list(edgelist$Gene1) %in% as.list(imap_comm_det[94])),]
#g<- graph.data.frame(imap_comm_det[1], directed=F)
#degree_centrality_module <- degree(g)


#predictive analysis per module
expression_data_gds5037 <- read.csv("exprsdata_gds5037.csv")
drops <- c("X")
expression_df_gds5037 <- expression_data_gds5037[ , !(names(expression_data_gds5037) %in% drops)]
metadata_gds5037 <- read.csv("metadata_gds5037.csv")
metadata_gds5037 <- metadata_gds5037[,-1]
require("randomForest")

rf_mod_i <- list()
rf_mod_j <- list()
rf_module_df <- list()
for(i in (1:total_com))
{
  module_i <- as.data.frame(imap_comm_det[i])
  if(nrow(module_i) > 10){
    print(i)
    colnames(module_i) <- 'ID_REF'
    for(j in 1:10){
      rf_module_df[[j]] <- module_i[sample(nrow(module_i),10),]
      rf_module_df[[j]] <- as.data.frame(rf_module_df[[j]])
      colnames(rf_module_df[[j]]) <- 'ID_REF'
      get_exp_data <- merge(expression_df_gds5037,rf_module_df[[j]], by = c('ID_REF'))
      drop_cols <- c("ID_REF", "IDENTIFIER")
      exp_df_for_rf <- get_exp_data[,!(names(get_exp_data) %in% drop_cols)]
      exp_df_for_rf_t <- t(exp_df_for_rf)
      exp_df_class <- cbind(metadata_gds5037$disease.state, exp_df_for_rf_t)
      exp_df_class <- as.data.frame(exp_df_class)
      exp_df_class$V1 <- as.factor(exp_df_class$V1)
      rf <- randomForest(V1~. , data=exp_df_class, ntree=1000,importance=T)
      print(rf$err.rate[1000])
      rf_mod_j[[j]] <- rf$err.rate[1000]

    }
    rf_mod_i[[i]] <- rf_mod_j
  }
}

rf_module_df <- as.data.frame(do.call(cbind, rf_mod_i))
write.csv(rf_module_df, "rf_modules.csv", row.names = F)
rf_csv <- read.csv("rf_modules.csv", header = T)
boxplot(rf_csv, xlab = "Modules", ylab = "Random forest error rate")
