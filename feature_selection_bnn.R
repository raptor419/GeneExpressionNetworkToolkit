#-------------------## feature selection BNN ## ------------------------#

#setwd("E:\\RShinyapp\\networksanalysis\\code")
#setwd("C://Networks//code//")
require(igraph)
require(BNN)
edgelist <- read.csv("Edgelist_top_10prcnt_genes.csv")
edgelist["abs_value"] <- abs(edgelist[,3])
edgelist = edgelist[,-3]
graph_obj <- graph.data.frame(edgelist)
imap_comm_det <- cluster_infomap(graph_obj, e.weights = edgelist[,3], v.weights = edgelist[,1], nb.trials = 10,
                                 modularity = F)
#total number of communities
total_com <- max(imap_comm_det$membership)
print(imap_comm_det$membership)
expression_data_gds5037 <- #edata
  read.csv("exprsdata.csv")
drops <- c("X")
expression_df_gds5037 <- expression_data_gds5037[ , !(names(expression_data_gds5037) %in% drops)]
diseased <- as.data.frame(read.csv("Diseased.csv"))
metadata_gds5037 <- #mdata
  read.csv("metadata.csv")
metadata_gds5037 <- metadata_gds5037[,-1]
require("randomForest")
require(Boruta)

bnn_features <- list()
bnn_names <- list()
for(i in (1:total_com))
{
  module_i <- as.data.frame(imap_comm_det[i])
  if(nrow(module_i) > 0){
    colnames(module_i) <- 'ID_REF'
    fetch_gene_expression<- merge(module_i,expression_df_gds5037, by='ID_REF',sort=F)
    rownames(fetch_gene_expression) <- fetch_gene_expression[,1]
    fetch_gene_expression <- fetch_gene_expression[,-(1:2)]
    rownames(fetch_gene_expression) <- as.factor(rownames(fetch_gene_expression))
    fetch_gene_expression_transpose <- t(fetch_gene_expression)
    fetch_gene_expression_class <- cbind(metadata_gds5037$disease.state, fetch_gene_expression_transpose)
    fetch_gene_expression_class <- as.data.frame(fetch_gene_expression_class)
    fg.X <- t(fetch_gene_expression[1:nrow(fetch_gene_expression),])
    ans <- colnames(fg.X)
    rownames(fg.X) <- NULL
    colnames(fg.X) <- NULL
    fg.Y <- diseased$V1
    bnn <- BNNsel(fg.X,fg.Y,total_iteration = 10000)
    for(j in (1:ncol(fg.X))){
      if(j<=length(bnn$fit) && bnn$fit[j]>0.5){
        bnn_features[[j]] <- ans[j]
        bnn_features[[j]] <- as.data.frame(bnn_features[[j]])
        print(ans[j])
      }
    }
    #png(paste("Boruta_module",i,".png"), width = 480, height = 480, units = "px", pointsize = 12)
    #plot(boruta,las=2, xlab= " ")
    #dev.off()
    #final decision tells confirmed rejected and tentative for all genes of a module
    #write.csv(boruta$finalDecision, paste("boruta_module",i,".csv"))
    
  }
}
require(plyr)
combine_bnn_features <- ldply(bnn_features,data.frame)
combine_bnn_features_omit_na <- na.omit(combine_bnn_features)
combine_bnn_features_omit_na <- as.data.frame(combine_bnn_features_omit_na)
if (length(bnn_features)!=0) {
  colnames(combine_bnn_features_omit_na) <- 'ID_REF'
  fetch_gene_expression_df<- merge(combine_bnn_features_omit_na,expression_df_gds5037, by='ID_REF',sort=F)
} else {
  combine_bnn_features_omit_na[1,1] <- 'No significant genes'
  colnames(combine_bnn_features_omit_na) <- 'ID_REF'
  fetch_gene_expression_df<- merge(combine_bnn_features_omit_na,expression_df_gds5037, by='ID_REF',sort=F)
}
write.csv(fetch_gene_expression_df,"Feature_Selection.csv")

