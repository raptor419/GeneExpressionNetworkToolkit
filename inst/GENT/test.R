

#norm_data()[,-((ncol(norm_data())-1): (ncol(norm_data())))]


expression_data_gds5037 <- #edata
  read.csv("exprsdata.csv")
drops <- c("X")
expression_df_gds5037 <- expression_data_gds5037[ , !(names(expression_data_gds5037) %in% drops)]
test <- expression_df_gds5037[,-(1:2)]
metadata_gds5037 <- #mdata
  read.csv("metadata.csv")
mdata <- metadata_gds5037
#rownames(test) <- edata[,2]
test[is.na(test)] <- as.numeric(0)
new<- t(test)
splitted_groups <- mdata[,1]
#give option to users to split by disease and cell type
required_format<- cbind(splitted_groups, new)
tempdf <- as.data.frame(required_format)
tempdf$splitted_groups <- as.factor(tempdf$splitted_groups)
splitted_groups <- split(tempdf,tempdf$splitted_groups)
if(input$diffexmethods == 'ano') {
  require(pbapply)
  pv <- pbapply(tempdf[,-1], 2, function(x){
    oneway.test(x ~ tempdf$splitted_groups,data=tempdf[,-1])$p.value
  })
  pvalue <- data.frame(pv)}
if(input$diffexmethods == 'eb'){
  m <- as.data.frame(tempdf[,-1])
  colnames(m) <- NULL
  rownames(m) <- NULL
  sapply(m, class)
  design <- as.numeric(tempdf[,1])
  fit <- lmFit(t(m),design = design)
  options(scipen = 999)
  pvalue <- eBayes(fit)$p.value
  pvalue <- data.frame(pvalue)}
if(input$diffexmethods == 'edr') {
  require(edgeR)
  f <- matrix(1, nrow = nrow(tempdf), ncol = 2)
  x <- 0:(nrow(tempdf)-1)
  f[,2]=x
  d <- DGEList(t(tempdf[,-1]))
  # Fit the NB GLMs
  disp <- estimateDisp(d, f)$common.dispersion
  fit <- glmQLFit(d ,f,dispersion = disp)
  # Likelihood ratio tests for trend
  pv <- glmLRT(fit)$table$PValue
  pvalue <- data.frame(pv)}
pvalue[is.na(pvalue)] <- as.numeric(1)
testdata <- cbind(pvalue,test)
colnames(testdata)[1] <- "pv"
deg<-subset(testdata,testdata$pv < 0.05)