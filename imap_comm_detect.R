## ------INFOMAP community detection-------
#setwd("E:\\RShinyapp\\networksanalysis\\code")
#setwd("C://Networks//code//")
require(igraph)
edgelist_el <- read.csv("Edgelist_top_10prcnt_genes.csv")
#edgelist <- edgelist_table()


# edgelist_el["abs_value"] <- abs(edgelist_el[,3])
# edgelist_el = edgelist_el[,-3]

graph_obj <- graph.data.frame(edgelist_el)
imap_comm_det <- cluster_infomap(graph_obj, e.weights = edgelist_el[,3], v.weights = edgelist_el[,1], nb.trials = 10,
                modularity = F)
#total number of communities
total_com <- max(imap_comm_det$membership)
#view each community
#imap_comm_det[1]
#imap_comm_det
