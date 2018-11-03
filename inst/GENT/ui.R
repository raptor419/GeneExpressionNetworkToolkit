#setwd("C://Networks//code//")
library(shiny)
library(shinyBS)
library(shinythemes)

shinyUI(navbarPage("GX Analyst", selected = "Show prioritised Gene List", theme = shinytheme("united"),
                   tags$style(".navbar{background-color:#B00000;}"),
                   tabPanel("Show prioritised Gene List",
                            
                            fluidRow(
                              column(4, wellPanel(theme = shinytheme("united"),
                                                  tags$style(".well {background-color:#B00000; color: #0A1432; border-color: #B00000;}"),
                                                  #tags$head(tags$style(".well{color: white;
                                                  #                          font-size: 15px;
                                                  #                   #font-style: italic;
                                                  #                  }")),
                                                  textInput("organism", "Enter organism", value = "Homo sapiens"),
                                                  textInput("disease", "Enter disease/cell line", value = "asthma"),
                                                  numericInput("geneids", "Minimum number of Gene-IDs per dataset", 20000),
                                                  #for coloring buttons use background-color, for text use color
                                                  #actionButton("goButton", "Get datasets", style = "background-color:#E8ADAA;"),
                                                  actionButton("goButton", "Get datasets"),
                                                  br(),
                                                  p("Displays all GEO public domain datasets with given specifications "),
                                                  br(),
                                                  textInput("gdsid", "Enter GDSID of experiment for downstream analysis", value = "GDS3615"),
                                                  actionButton("expbutton", "View expression data"),
                                                  actionButton("metbutton", "View metadata"),
                                                  downloadLink('downloadData', 'Download data'),
                                                  hr(style = "border-top: 1px solid white;"),
                                                  h4(p("Data Normalization")),
                                                  selectInput("normmethods", "Choose normalization method:",
                                                              choices = c('Log'='log','Quantile'='qtl', 'Log+Quantile' = 'lqt', 'Z-score' = 'zs', 'Unitization' = 'lw', 'Normalisation about median' = 'qsp'), selected = 'lqt'),
                                                  actionButton("normbutton", "Normalize data"),
                                                  actionButton("normplotbutton", "View plots"),
                                                  hr(),
                                                  h4(p("Minimal Gene List")),
                                                  actionButton("minimal", "Show genes")
                              )),
                              
                              column(8,tabsetPanel
                                     (
                                       #tags$style(".{background-color:#7F525D; color: #EDC9AF; border-color: #7D0552;}"),
                                       tabPanel('View All datasets',dataTableOutput("download_df")),
                                       tabPanel('View Expression data',dataTableOutput("expdata")),
                                       tabPanel('View Metadata',dataTableOutput("metdata")),
                                       tabPanel('View Normalized data', dataTableOutput("norm_df")),
                                       tabPanel('Show prioritised gene list',dataTableOutput("mini_gl")),
                                       
                                       mainPanel(
                                         bsModal("modalExample", "Raw vs Normalized Data Plots", "normplotbutton", size = "large",plotOutput("plot1"), plotOutput("plot2"),downloadButton('downloadPlot', 'Download'))
                                       )
                                       
                                     ))
                              
                              
                            )
                   ),
                   ##---------------------------------CUSTOM ANALYSIS TAB -------------------------------###
                   tabPanel("Customise Process",
                            fluidRow(
                              column(4, wellPanel(
                                h4("Differential Gene Expression"),
                                #h5("Choose dataset"),
                                checkboxInput("raw_df", label = "Use raw data", value = FALSE),
                                checkboxInput("norm_df", label = "Use normalized data", value = TRUE),
                                textInput("pval", "P value <=", value = "0.05"),
                                checkboxInput("checkbox_fdr", label = "False discovery rate (fdr)", value = TRUE),
                                selectInput("diffexmethods", "Choose differential expression calculation method:",
                                            choices = c('ANOVA'='ano','eBayes(Recommended)'='eb', 'edgeR'='edr'),selected = 'eb'),
                                actionButton("diffexpbutton", "Get differentially expressed genes"),
                                h4("Association Analysis"),
                                selectInput("assindex", "Choose association index:",
                                            choices = c('Pearson Correlation'='pcor',
                                                        'Maximal Information Nonparametric Exploration(MINE)(Recommended)'='mine',
                                                        'Jaccard Index' = 'jac', 'Cosine' = 'cos', 'Sorensen' = 'simp',
                                                        'Minkowski' = 'geo','Spearman' = 'hgeo'), selected = 'pcor'),
                                actionButton("assoc_but", "Get Association matrix"),
                                textInput("frac_genes", "Pick top X-percent gene relationships (Enter X)", value = "10"),
                                p("Based on association index, an egdelist is created, out which top 10% genes with
                                  highest association index are selected for downstream analysis"),
                                actionButton("edglist_but", "Get Edgelist"),
                                hr(),
                                h4("Gene Network"),
                                actionButton("network_but", "View gene network"),
                                h4("Using InfoMAP community detection algorithm"),
                                actionButton("net_chars_but", "Get Module Characteristics"),
                                textInput("module_num", "Select module number to view genes in it", value = "1"),
                                actionButton("map_but", "Get module"),
                                hr(),
                                h4("Predictive Analysis"),
                                strong("Network Topology Analysis"),
                                actionButton("mod_rf_but", "Plot-Predictive power of modules"),
                                p("Find out significant modules in network based on its predictive accuracy"),
                                #strong("Using degree centraliity as network parameter"),
                                #actionButton("plot_degcen_but", "Plot"),
                                # p("Plots predictive power vs degree centrality"),
                                h4("Feature selected genes using Bayesian Neural Networks"),
                                #selectInput("fs", "Feature selected genes using Bayesian Neural Networks",
                                #choices = c("Boruta" = "bor","Bayesian Neural Networks" = "bnn"), selected = "bor"),
                                actionButton("fs_but", "Feature selection output"),
                                hr(),
                                h4("Final Gene List"),
                                strong("Minimal set of differentially expressed genes with discriminatory predictive power"),
                                actionButton("gene_list_but", "View gene list"),
                                br(),
                                strong("Compare predictive power of selected gene list vs diff exp gene list"),
                                selectInput("ppc", "Choose algorithm for comparison",
                                            choices = c("Boruta" = "bor","Bayesian Neural Networks" = "bnn"), selected = "bor"),
                                actionButton("plot_fs","Predictive power comparison plot" ),
                                p("Boxplot comparing the predictive power of feature selected genes vs top differentially expressed genes")
                                #closing bracket for column4
                                )),
                              column(8,tabsetPanel(
                                tabPanel('View differentially expressed genes',dataTableOutput("diff_exp")),
                                tabPanel('View Association Matrix',dataTableOutput("ass_matrix")),
                                tabPanel('View Edgelist',dataTableOutput("edgelist")),
                                tabPanel('View Community characteristics',dataTableOutput("cec")),
                                tabPanel("View Modules", dataTableOutput("map")),
                                tabPanel("View Feature Selected Genes", dataTableOutput("fs_bor_df")),
                                tabPanel("View Final Gene List", dataTableOutput("final_genes")),
                                mainPanel(theme = shinytheme("united"),
                                          bsModal("modalExample2", "Gene Network", "network_but", size = "large",
                                                  htmlOutput("network_graph"), downloadButton('downloadnetwork', 'Save HTML') ),
                                bsModal("modalExample3", "Random Forest error rate between gene communities",
                                        "mod_rf_but", size = "large",plotOutput("rf_module_plot", "100%", "500px"),downloadButton('download_rf_Plot', 'Save Plot')),
                                #bsModal("modalExample4", "Random Forest error rate wrt network parameters",
                                # "plot_degcen_but", size = "large",plotOutput("rf_degcen_plot", "100%", "500px"),downloadButton('download_rf_degcen', 'Save Plot')),
                                bsModal("modalExample5", "Predictive power plot comparison",
                                        "plot_fs", size = "large",imageOutput("fs_comp_plot", "100%", "500px"),downloadButton('download_comparison', 'Save Plot'))
                              )
                              ))
                              
                              
                              #closing bracket for fluid row
                            )
                            #clossing bracket of networks tab
                   ),
                   #tabPanel("Transcription factor analysis"),
                   tabPanel("About GX Analyst"),
                   tabPanel("About Us")
                   
))
