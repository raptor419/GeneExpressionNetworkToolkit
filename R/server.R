#setwd("C://Networks//code//")
# Make sure GEOmetadb.sqlite is present in the working directory. Also update the GEOsqlite 
# file to the latest version
# other libraries install pbapply, limma, shinyBS,shinythemes, shiny, networkd3,
# change working directory in all code files
## For first time running: 
  # Install libraries 
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("BiocUpgrade")
  # biocLite("GEOmetadb")
  # biocLite("Biobase")
  # biocLite("limma")
  # biocLite("edgeR")
  # getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz")
  # install.packages("ggplot2")
  # install.packages("philentropy")
  # install.packages("minerva")
  # install.packages("BNN")
  # install.packages("data.table")
  # install.packages("plotly")
  # install.packages("pbapply")
  # install.packages("reshape2")
  # install.packages("igraph")


library(GEOmetadb)
library(Biobase)


server <- function(input, output) {
  # Takes an action every time button is pressed;
  # printing a message to the console for log history
  observeEvent(input$goButton, {
    cat("Showing max genes", input$geneids, "rows\n")
  })
  # Take a reactive dependency on input$button, but
  # not on any of the stuff inside the function
  
  library(GEOmetadb)
  library(Biobase)
  library(ggplot2)
  library(philentropy)
  library(minerva)
  library(MEGENA)
  library(BNN)
  
  #---downlaod the latest version of GEOdb---#
  #getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz")
  
  con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
  query_str1 <- paste("SELECT * FROM gds WHERE sample_organism LIKE '%")
  query_str2 <- paste("%' AND sample_type LIKE '%RNA")
  query_str3 <- paste("%' AND description LIKE '%")
  query_str4 <- paste("%'")
  
  # view all datasets (GDSIDs) for input keyword
  df <- eventReactive(input$goButton, {
    dio <- paste(input$organism)
    dis <- paste(input$disease)
    query <- paste(query_str1, dio, query_str2, query_str3, dis, query_str4,sep="")
    gds_subset <- dbGetQuery(con, query)
    g <- subset(gds_subset,gds_subset$feature_count >= input$geneids)
    library(data.table)
    a <- c("timecourse","time series","time-series","timeseries","timepoint","time point","time-point","time-course","time course")
    g[- grep(paste(a,collapse="|"),g[,3]), ] 
  })
  output$download_df <- renderDataTable({df()})
  
  #view expression data
  expressiondata <- eventReactive(input$expbutton, {
    #options('download.file.method.GEOquery'='auto')
    dio <- paste(input$organism)
    dis <- paste(input$disease)
    query <- paste(query_str1, dio, query_str2, query_str3, dis, query_str4,sep="")
    print(query)
    gds_subset <- dbGetQuery(con, query)
    gds_subset <- subset(gds_subset,gds_subset$feature_count >= input$geneids)
    if (interactive()) {  
      #download data --  //// 
      output$downloadData <- downloadHandler(
        filename = function() {
          paste('data-', Sys.Date(), '.csv', sep='')
        },
        content = function(con) {
          gds_id<-gds_subset[,2]
          gds <- getGEO(gds_id[which(gds_id == input$gdsid)])
          mdata <- Columns(gds)
          edata <- Table(gds)
          #write.csv(mdata, con)
          write.csv(edata, con)
        })
    }
    gds_id<-gds_subset[,2]
    gds <- getGEO(gds_id[which(gds_id == input$gdsid)])
    Table(gds)
    
  })
  output$expdata <- renderDataTable({expressiondata()})
  
  #view metadata
  metadata <- eventReactive(input$metbutton, {
    dio <- paste(input$organism)
    dis <- paste(input$disease)
    query <- paste(query_str1, dio, query_str2, query_str3, dis, query_str4,sep="")
    print(query)
    gds_subset <- dbGetQuery(con, query)
    gds_subset <- subset(gds_subset,gds_subset$feature_count >= input$geneids)
    gds_id<-gds_subset[,2]
    gds <- getGEO(gds_id[which(gds_id == input$gdsid)])
    Columns(gds)})
  output$metdata <- renderDataTable({metadata()})
  
  #normalization
  norm_data <- eventReactive(input$normbutton, { 
    dio <- paste(input$organism)
    dis <- paste(input$disease)
    query <- paste(query_str1, dio, query_str2, query_str3, dis, query_str4,sep="")
    print(query)
    gds_subset <- dbGetQuery(con, query)
    gds_subset <- subset(gds_subset,gds_subset$feature_count >= input$geneids)
    gds_id<-gds_subset[,2]
    gds <- getGEO(gds_id[which(gds_id == input$gdsid)])
    mdata <- Columns(gds)
    edata <- Table(gds)
    resul <- NULL
    params <- NULL
    if(input$normmethods == 'log') {
      transform_data<- log(edata[,(3:ncol(edata))],base=2)}
    if(input$normmethods == 'qtl'){
      library(limma)
      transform_data <- normalizeQuantiles(edata[,3:ncol(edata)])}
    if(input$normmethods == 'zs'){
      transform_data <- as.data.frame(scale(edata[,3:ncol(edata)]))}
    if(input$normmethods == 'lw'){
      x <- edata[,3:ncol(edata)]
      for (i in 1:ncol(x)) {
        resul <- cbind(resul, (x[, i] - mean(x[, i]))/(max(x[, i]) - min(x[, i])))}
      transform_data <- as.data.frame(resul)}
    if(input$normmethods == 'qsp'){
      x <- edata[,3:ncol(edata)]
      for (i in 1:ncol(x)) {
        resul <- cbind(resul, (x[, i] - median(x[, i]))/(mad(x[, i])))}
      transform_data <- as.data.frame(resul)}
    if(input$normmethods == 'lqt'){
      log_transform<- log(edata[,(3:ncol(edata))],base=2)
      library(limma)
      transform_data <- normalizeQuantiles(log_transform)}
    transform_data["ID_REF"] = edata$ID_REF
    transform_data["IDENTIFIER"] = edata$IDENTIFIER
    rownames(transform_data) <- transform_data$ID_REF
    a <- expressiondata()[,-(1:2)]
    transform_data
  })
  output$norm_df <- renderDataTable({norm_data()})  
  
  library(plotly)
  #View plots
  output$plot1 <- 
    #renderPlot({ggplotly(norm_data()[,-((ncol(norm_data())-1): (ncol(norm_data())))], aes(x=cond, y=rating))+geom_boxplot()},height = 400,width = 600)
    renderPlot({boxplot(expressiondata()[,-(1:2)], main = "Raw expression data")})
  
  output$plot2 <- renderPlot({boxplot(norm_data()[,-((ncol(norm_data())-1): (ncol(norm_data())))], main = "Normalized data")})
  
  # download button for raw vs norm plots ////
  plotInput <- function(){boxplot(norm_data())}
  output$downloadPlot <- downloadHandler(
    filename = function() { paste('Shiny', '.pdf', sep='') },
    content = function(file) {
      pdf(file)
      boxplot(expressiondata()[,-(1:2)], main = "Raw expression data")
      boxplot(norm_data()[,-((ncol(norm_data())-1): (ncol(norm_data())))], main = "Normalized data")
      #print(plotInput())
      #plotInput2()
      dev.off()
      #ggsave(png(file), plot = plotInput())
    }) 
  
  #Minimal gene List
  min_gl <- eventReactive(input$minimal, { 
    e <- new.env()
    e$edata <- norm_data()
    e$idf <- input$gdsid
    #read.csv("exprsdata_gds5037.csv")  ##expressiondata()
    e$mdata <- metadata()
    diseased <- matrix(metadata()[,2])
    s <- as.data.frame(levels(metadata()[,2]))
    for(i in (1:nrow(metadata())))
    {
      diseased[i] <- which(grepl(metadata()[i,2], s$`levels(metadata()[, 2])`)) - 1
    }
    e$diseased <- diseased
    #read.csv("metadata_gds5037.csv") ##metadata()
    sys.source("minimal_gl.R", e)
    #feature_selection_output()
    e$fetch_gene_expression_df
  })
  
  
  output$mini_gl <- renderDataTable({min_gl()})
  
  ##-----------------------CUSTOM ANALYSIS TAB------------------------------------##
  
  
  
  #differentially expressed genes
  diffexp_df <- eventReactive(input$diffexpbutton, {
    edata <- expressiondata()
    mdata <- metadata()
    if(input$raw_df == TRUE){ test <- edata[,-(1:2)]}
    if(input$norm_df == TRUE){test <- norm_data()[,-((ncol(norm_data())-1): (ncol(norm_data())))]}

    #rownames(test) <- edata[,2]
    #test <- test[(1:500),]
    test[is.na(test)] <- as.numeric(0)
    new<- t(test)
    splitted_groups <- mdata[,2]
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
      pvalue <- data.frame(as.numeric(pv))}
    pvalue[is.na(pvalue)] <- as.numeric(1)
    
    #colnames(pvalue)
    #false discovery rate:-
    if(input$checkbox_fdr==TRUE){
      fdr <- p.adjust(pvalue[,1], method="fdr")
      testdata <- cbind(fdr,test)
      deg<-subset(testdata,testdata$fdr < input$pval)
    }
    else{
      testdata <- cbind(pvalue,test)
      colnames(testdata)[1] <- "pv"
      deg<-subset(testdata,testdata$pv < input$pval)
      #deg[,-1]
    }
    #deg[,-1]
    deg["ID_REF"] = rownames(deg)
    deg
  })
  output$diff_exp <- renderDataTable({diffexp_df()})  
  
  ##get diff exp dataset
  ## get association matrix
  #Edgelist
  assoc_matrix <- eventReactive(input$assoc_but, {
    deg <- diffexp_df()
    drop <- c("fdr", "ID_REF")
    deg_droped = deg[,!(names(deg) %in% drop)]
    if(input$assindex == 'pcor') {
      assoc <- cor(t(deg_droped), method = "pearson")
    }
    if(input$assindex == 'mine') {
      assoc <- mine(t(deg_droped))$MIC
      colnames(assoc) <- colnames(t(deg_droped))
    }
    if(input$assindex == 'jac') {
      assoc <- distance(deg_droped, method = "jaccard")
      colnames(assoc) <- colnames(t(deg_droped))
    }
    if(input$assindex == 'cos') {
      assoc <- distance(deg_droped, method = "cosine")
      colnames(assoc) <- colnames(t(deg_droped))
    }
    if(input$assindex == 'simp') {
      assoc <- distance(deg_droped, method = "sorensen")
      colnames(assoc) <- colnames(t(deg_droped))
    }
    if(input$assindex == 'geo') {
      assoc <- distance(deg_droped, method = "minkowski")
      colnames(assoc) <- colnames(t(deg_droped))
    }
    if(input$assindex == 'hgeo') {
      assoc <- cor(t(deg_droped), method = "spearman")
    }
    assoc})
  output$ass_matrix <- renderDataTable({assoc_matrix()})
  
  #correlation
  #wgcna
  #adjacency = adjacency(datExpr, power = softPower)
  #-----------------------------------------------------------------------------------------#
  #Edgelist
  edgelist_table <- eventReactive(input$edglist_but, {
    #mi_nonzero = adjacency(assoc_matrix(), power = softPower)
    ab <- input$frac_genes
    mi_nonzero <- assoc_matrix()
    diag(mi_nonzero) <- 0
    require(reshape2)
    edgelist_mi_nonzero <- melt(mi_nonzero)
    edglist <- calculate.PFN(as.data.frame(edgelist_mi_nonzero))
    ac <- nrow(edglist)
    ad <- as.numeric(ab)*as.numeric(ac)/100
    edglist <- edglist[1:ad,]
    #require(reshape2)
    #edgelist_mi_nonzero <- melt(mi_nonzero)
    #dec_order <- edgelist_mi_nonzero[order(edgelist_mi_nonzero$value , decreasing = TRUE),]
    #dec_order <- subset(dec_order,dec_order$value!=0) 
    #edglist <- calculate.correlation(dec_order,method = "pearson",FDR.cutoff = 0.05)
    #per = as.numeric(input$frac_genes)*nrow(dec_order))
    #el <- dec_order[(1:((as.numeric(input$frac_genes))/100*(nrow(dec_order)))),]
  })
  output$edgelist <- renderDataTable({edgelist_table()})
  
  #-------------------------------------------------------------------------------# 
  saveData <- function(data) {
    diseased <- matrix(metadata()[,2])
    s <- as.data.frame(levels(metadata()[,2]))
    for(i in (1:nrow(metadata())))
    {
      diseased[i] <- which(grepl(metadata()[i,2], s$`levels(metadata()[, 2])`)) - 1
    }
    write.csv(diseased,row.names = FALSE, "Diseased.csv")
    write.csv(edgelist_table(),row.names = FALSE, "Edgelist_top_10prcnt_genes.csv")
    write.csv(metadata(),row.names = FALSE, "metadata.csv")
    write.csv(expressiondata(),row.names = FALSE, "exprsdata.csv")
    write.csv(diffexp_df(),row.names = FALSE, "Diff_Exp_Data.csv")
  }
  
  #Gene Network Plot
  getNetwork<-function() {
    saveData(edgelist_table)
    return(source("networkd3.R"))
  }
  output$network_graph<-renderUI({getNetwork()})
  
  #------------------------------------------------------------------------------#
  #View modules
  
  #infomap <- function(x){
  #return(source("imap_comm_detect.R"))
  # return(source("imap_comm_detect.R", local = TRUE))
  #}
  comm_det_modules <- eventReactive(input$map_but, {
    require(igraph)
    e <- new.env() ## new.env so exchange variables btw source file and main file
    e$edgelist_el <- edgelist_table()
    sys.source("imap_comm_detect.R", e)
    #infomap(edgelist_el)
    module_view <- as.data.frame(e$imap_comm_det[input$module_num])
    module_view
  })
  output$map <- renderDataTable({comm_det_modules()})  
  #-------------------------------------------------------------------------#
  
  
  #Network characteristics
  net_chars <- eventReactive(input$net_chars_but, {
    require(igraph)
    e <- new.env() ## new.env so exchange variables btw source file and main file
    e$edgelist_el <- edgelist_table()
    sys.source("imap_characteristics.R", e)
    #infomap(edgelist_el)
    aj <- as.data.frame(e$comm_list)
    aj
  })
  output$cec <- renderDataTable({net_chars()})
  #------------------------------------------------------------------------#  
  ##RF error rate between modules
  #rf_modules <- function(){return(source("random_forest_modules.R"))}
  
  rf_modules_df <- eventReactive(input$mod_rf_but, {
    e <- new.env()
    e$edata <- norm_data()
    #read.csv("exprsdata.csv")  ##expressiondata()
    e$mdata <- metadata()
    #read.csv("metadata.csv") ##metadata()
    e$edgelist <- edgelist_table()
    sys.source("random_forest_modules.R", e)
    e$rf_csv
    #rf_modules()
    #rf_csv
  })
  
  output$rf_module_plot <- renderPlot({
    input$mod_rf_but
    boxplot(rf_modules_df(), xlab = "Modules", ylab = "Random forest error rate")
  })
  
  # download button for RF b/w modules plot- ////
  plotInput <- function(){
    boxplot(rf_modules_df())}
  output$download_rf_Plot <- downloadHandler(
    filename = function() { paste('RF_btw_modules', '.pdf', sep='') },
    content = function(file) {
      png(file)
      boxplot(rf_modules_df())
      dev.off()}) 
  #------------------------------------------------------------------------#
  
  #RF error rate vs network parameters
  
  #degree centrality
  
  #rf_degree_centrality <- function(){
  #return(source("random_forest_network_params.R"))
  #}
  rf_deg_cent_res <- eventReactive(input$plot_degcen_but, {
    e <- new.env()
    e$edata <- norm_data()
    #read.csv("exprsdata_gds5037.csv")  ##expressiondata()
    e$mdata <- metadata()
    #read.csv("metadata_gds5037.csv") ##metadata()
    e$edgelist <- edgelist_table()
    if(input$netpars == 'dcen'){
      e$net_param = 'degree_cent'
    }
    else if(input$netpars == 'bet'){
      e$net_param = "edge_bet"
    }
    #e <- new.env()
    #e$edata <- read.csv("exprsdata_gds5037.csv")  ##expressiondata()
    #e$mdata <- read.csv("metadata_gds5037.csv") ##metadata()
    #e$edgelist <- edgelist_table()
    sys.source("random_forest_network_params.R", e)
    e$rf_result
    #rf_degree_centrality()
    #rf_result
    
  })
  
  output$rf_degcen_plot <- renderPlot({
    input$plot_degcen_but
    rf_result_df <- as.data.frame(rf_deg_cent_res())
    rf_result_reqd <- cbind(rownames(rf_result_df),rf_result_df)
    rf_result_numeric <- apply(rf_result_reqd, 2 , as.numeric)
    m <- lm(rf_result_numeric[,2] ~ rf_result_numeric[,1])
    plot(rf_deg_cent_res(), type="p",main= bquote("Slope"== .(m$coefficients[2])),  
         xlab="Models with probes of high degree to low", ylab="Error rate of RF model")
    abline(m, col="red")
  })
  
  # download button for RF degree centrality plot- ////
  plotInput <- function(){
    plot(rf_deg_cent_res(), type="p",main= "Out Of Bag Error",  
         xlab="Models with probes of high degree to low", ylab="Error rate of RF model")
  }
  output$download_rf_degcen <- downloadHandler(
    filename = "RF_deg_cent.png",
    content = function(file) {
      png(file)
      plot(rf_deg_cent_res(), type="p",main= "Out Of Bag Error",  
           xlab="Models with probes of high degree to low", ylab="Error rate of RF model")
      dev.off()}) 
  
  #-----------------------------------------------------------------------------------#
  
  #Feature selection 
  
  #feature_selection_output <- function(){
  #return(source("feature_selection.R"))
  #}
  
  #feature selection output df with conf, tent and rejected gene exp data info
  fs_genes_df <- eventReactive(input$fs_but, {
    e <- new.env()
    e$edata <- norm_data()
    e$idf <- input$gdsid
    #read.csv("exprsdata_gds5037.csv")  ##expressiondata()
    e$mdata <- metadata()
    #read.csv("metadata_gds5037.csv") ##metadata()
    e$edgelist <- edgelist_table()
      sys.source("feature_selection_bnn.R", e)
    #feature_selection_output()
    e$fetch_gene_expression_df
  })
  output$fs_bor_df <- renderDataTable({fs_genes_df()})
  #-----------------------------------------------------------------------------------#
  #Final confirmed gene list(with gene symbols) from feature selection bor
  #taken selected genes from each module and combined
  final_gene_df <- eventReactive(input$gene_list_but, {
    fs_genes_df()[,(1:2)]
    #gene_df
  })
  output$final_genes <- renderDataTable({final_gene_df()})
  
  #------------------------------------------------------------------------------------#
  #Boxplot for comparison of predictive accuracy of selected genes 
  #vs differentially expressed genes
  
  rf_modules_d <- eventReactive(input$plot_fs, {
    e <- new.env()
    e$edata <- norm_data()
    #read.csv("exprsdata.csv")  ##expressiondata()
    e$mdata <- metadata()
    #read.csv("metadata.csv") ##metadata()
    e$edgelist <- edgelist_table()
    if(input$ppc == 'bor') {
      e$algo <- 'bor'
    }
    if(input$ppc == 'bnn') {
      e$algo <- 'bnn'
    }
    sys.source("random_forest_genes.R", e)
    e$rf_csv
  })
  
  output$fs_comp_plot <- renderPlot({
    input$plot_fs
    boxplot(rf_modules_d(), xlab = "Genes",names = c("Differentially expressed genes" ,"Feature selected genes"), ylab = "Random forest error rate")
  })
  
  # download button for RF b/w modules plot- ////
  plotInput <- function(){
    boxplot(rf_modules_d())}
  output$download_comparison <- downloadHandler(
    filename = "RF_Genes.png",
    content = function(file) {
      png(file)
      #plotInput()
      boxplot(rf_modules_d())
      dev.off()}) 
  
  
  
  
  #last closing bracket    
}

