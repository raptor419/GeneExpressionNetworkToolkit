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
net_param = "degree_cent"
expression_data_gds5037 <- read.csv("exprsdata.csv")
drops <- c("X")
expression_df_gds5037 <- expression_data_gds5037[ , !(names(expression_data_gds5037) %in% drops)]
metadata_gds5037 <- read.csv("metadata.csv")
metadata_gds5037 <- metadata_gds5037[,-1]
diseased <- as.data.frame(read.csv("Diseased.csv"))

require("randomForest")
#degree centrality for genes in all modules window sliding down top to bottom overlap--
degree_centrality_module_dec <- list()
Edgelist_cl <- list()


for(i in (1:total_com))
{
  if(nrow(as.data.frame(imap_comm_det[i])) > 0){
    cond1<-edgelist[,1] %in% as.data.frame(imap_comm_det[i])[,1]
    cond2<-edgelist[,2] %in% as.data.frame(imap_comm_det[i])[,1]

    cond3<-1:length(cond1)
    for(j in 1:length(cond3)){
      if(cond1[j]==TRUE && cond2[j]==TRUE){
        cond3[j]<-TRUE
      }
      else{
        cond3[j]<-FALSE
      }
    }
    Edgelist_cl[[i]] <- edgelist[cond3==T,]
    g<- graph.data.frame(Edgelist_cl[[i]], directed=F)

    if(net_param == "degree_cent"){
      degree_centrality_module <- degree(g)}
    else if(net_param == "edge_bet"){
      degree_centrality_module <- degree(g)
    }

    degree_centrality_module <- as.data.frame(degree_centrality_module)
    degree_centrality_module_genes <- cbind(rownames(degree_centrality_module),degree_centrality_module)

    degree_centrality_module_dec[[i]] <- degree_centrality_module_genes[order(degree_centrality_module_genes$degree_centrality_module,decreasing=T),]
    degree_centrality_module_dec[[i]] <- as.data.frame(degree_centrality_module_dec[[i]])
    }
}

###-----random forest of modules-----##
test_data_rf <- degree_centrality_module_dec
print(test_data_rf)

require(randomForest)
rf_result <- rep(0,10)
for(i in 1:10){
  genes <- c()

  for(j in 1:total_com) #length(test_data_rf)
  {
    if(nrow(test_data_rf[[j]]) > 0){
      length_module_to_be_taken <- nrow(test_data_rf[[j]])
      whole_number_module_genes <- test_data_rf[[j]][(1:length_module_to_be_taken),1]
      part_length <- length(whole_number_module_genes)/10
      g <- factor(rep(1:10, each=part_length))
      genes <- c(genes,as.vector(split(whole_number_module_genes,g)[i]))

    }
  }
  genes <- data.frame(unlist(genes))
  #write.csv(genes,paste("subset",i,".csv"))
  colnames(genes) <- 'ID_REF'
  rownames(genes) <- genes[,1]
  fetch_gene_expression<- merge(genes,expression_df_gds5037, by='ID_REF',sort=F)
  fetch_gene_expression <- fetch_gene_expression[,-(1:2)]

  fetch_gene_expression_transpose <- t(fetch_gene_expression)
  #for all datasets= mydata[,3] , for copd = mydata[,2]
  fetch_gene_expression_class <- cbind(metadata_gds5037$disease.state, fetch_gene_expression_transpose)
  fetch_gene_expression_class <- as.data.frame(fetch_gene_expression_class)
  rf <- randomForest(diseased$x~. , data=fetch_gene_expression_class, ntree=1000,importance=T)
  #rf_result[i] <- rf$test
  rf_result[i]<- rf$err.rate[1000]
}
#plot(rf_result, type="l",main= "Out Of Bag Error", xlab="Models with probes of high degree to low", ylab="Error rate of RF model")
#plot(rf_result, type="p",main= "Out Of Bag Error",  xlab="Models with probes of high degree to low", ylab="Error rate of RF model")
#rf_result_df <- as.data.frame(rf_result)
#rf_result_reqd <- cbind(rownames(rf_result_df),rf_result_df)
#rf_result_numeric <- apply(rf_result_reqd, 2 , as.numeric)
#m <- lm(rf_result_numeric[,2] ~ rf_result_numeric[,1])
#abline(m, col="red")
