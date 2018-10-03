-------------------## feature selection Boruta ## ------------------------

#setwd("E:\\RShinyapp\\networksanalysis\\code")
#setwd("C://Networks//code//")
require(igraph)
edgelist <- read.csv("Edgelist_top_10prcnt_genes.csv")
graph_obj <- graph.data.frame(edgelist)
imap_comm_det <- cluster_infomap(graph_obj, e.weights = edgelist[,3], v.weights = edgelist[,1], nb.trials = 10,
                                 modularity = TRUE)
#total number of communities
total_com <- max(imap_comm_det$membership)

expression_data_gds5037 <- read.csv("exprsdata.csv")
drops <- c("X")
expression_df_gds5037 <- expression_data_gds5037[ , !(names(expression_data_gds5037) %in% drops)]
metadata_gds5037 <- read.csv("metadata.csv")
metadata_gds5037 <- metadata_gds5037[,-1]
require("randomForest")



require(Boruta)
boruta_features <- list()
for(i in (1:total_com))
{
  module_i <- as.data.frame(imap_comm_det[i])
  if(nrow(module_i) > 40){
    colnames(module_i) <- 'ID_REF'
    fetch_gene_expression<- merge(module_i,expression_df_gds5037, by='ID_REF',sort=F)
    rownames(fetch_gene_expression) <- fetch_gene_expression[,1]
    fetch_gene_expression <- fetch_gene_expression[,-(1:2)]
    rownames(fetch_gene_expression) <- as.factor(rownames(fetch_gene_expression))
    fetch_gene_expression_transpose <- t(fetch_gene_expression)
    fetch_gene_expression_class <- cbind(metadata_gds5037$disease.state, fetch_gene_expression_transpose)
    fetch_gene_expression_class <- as.data.frame(fetch_gene_expression_class)
    fetch_gene_expression_class$V1 <- as.factor(fetch_gene_expression_class$V1)
    boruta <- Boruta(V1~. , data=fetch_gene_expression_class)
    #png(paste("Boruta_module",i,".png"), width = 480, height = 480, units = "px", pointsize = 12)
    plot(boruta,las=2, xlab= " ")
    #dev.off()
    boruta_features[[i]] <- getSelectedAttributes(boruta)
    #boruta_features[[i]] <- as.data.frame(boruta_features[[i]])
    #final decision tells confirmed rejected and tentative for all genes of a module
    #write.csv(boruta$finalDecision, paste("boruta_module",i,".csv"))

}
}


#Borute confirmed genes dataframe with expression values, identifier, etc
require(plyr)
combine_boruta_features <- ldply(boruta_features,data.frame)
combine_boruta_features_omit_na <- na.omit(as.vector(unlist(combine_boruta_features)))
combine_boruta_features_omit_na <- as.data.frame(combine_boruta_features_omit_na)
combine_boruta_features_omit_na[,1] <- gsub("`", '', combine_boruta_features_omit_na[,1])
rownames(combine_boruta_features_omit_na) <- as.factor(combine_boruta_features_omit_na[,1])
colnames(combine_boruta_features_omit_na) <- 'ID_REF'
fetch_gene_expression_df<- merge(combine_boruta_features_omit_na,expression_df_gds5037, by='ID_REF',sort=F)
#write.csv(fetch_gene_expression,"Boruta_gene_expression.csv")


#Final gene list
gene_df <- fetch_gene_expression_df[,(1:2)]

#Rf error rate of borute confirmed genes for comparison plot
fetch_gene_expression_rf <- fetch_gene_expression_df[,-(1:2)]
fetch_gene_expression_transpose <- t(fetch_gene_expression_rf)
fetch_gene_expression_class <- cbind(metadata_gds5037$disease.state, fetch_gene_expression_transpose)
fetch_gene_expression_class <- as.data.frame(fetch_gene_expression_class)
fetch_gene_expression_class$V1 <- as.factor(fetch_gene_expression_class$V1)
rf <- randomForest(V1~. , data=fetch_gene_expression_class, ntree=10000,importance=T)
rf_result_boruta_genes<- rf$err.rate[10000]

#Rf error rate of differentially expressed genes

# plot rf_result_boruta_genes, rf_results_diff_genes
#write code here
