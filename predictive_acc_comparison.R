#Predictive accuracy comparison plot
#setwd("C://Networks//code//")

#fetch_gene_expression_df can be taken from feature selection script
metadata_gds5037 <- mdata

#Rf error rate of borute confirmed genes for comparison plot
fetch_gene_expression_rf <- fetch_gene_expression_df[,-(1:2)]
fetch_gene_expression_transpose <- t(fetch_gene_expression_rf)
fetch_gene_expression_class <- cbind(metadata_gds5037$disease.state, fetch_gene_expression_transpose)
fetch_gene_expression_class <- as.data.frame(fetch_gene_expression_class)
fetch_gene_expression_class$V1 <- as.factor(fetch_gene_expression_class$V1)
rf <- randomForest(V1~. , data=fetch_gene_expression_class, ntree=10000,importance=T)
rf_result_boruta_genes<- rf$err.rate[10000]

#Rf error rate of differentially expressed genes

#diff_exp_genes_df <-
rf_result_randomly_picked_genes <- c(0.67,0.61,0.65,0.55)
#see scale_up_module_detection

# plot rf_result_boruta_genes, rf_results_diff_genes
boxplot(rf_result_randomly_picked_genes,ylab= "Random Forest Error rate", ylim = c(0.2,0.8))
points(0.37,pch=19,col="red",cex=1)
