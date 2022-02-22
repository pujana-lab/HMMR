library(dplyr)
library(stringr)
library(network)
library(ggnetwork)
library(ggplot2)

# This script computes a weighted score for each bidirectional pathway and plots it.

# Read the file (a boolean matrix of keywords for each pathway)
keywords.df <- read.table("./Pascal_genesets_keywords.txt",
                          sep="\t", stringsAsFactors = F, header = T, row.names = 1)

# Get all keywords per geneset
keywords.df <- apply(keywords.df, 1, FUN=function(x){colnames(keywords.df)[x]})

# Remove genesets with no keywords
keywords.df <- keywords.df[lapply(keywords.df,length)>0]

# Get all combinations per geneset and their individual weight
keywords.w <- lapply(keywords.df, FUN=function(x){expand.grid(x,x,weight=1/length(x))})

# Collapse all genesets and remove pairs with same keywords
keywords.w <- bind_rows(keywords.w) %>% dplyr::filter(Var1 != Var2)

# Change the bidirectional pairs to a same order (to sum later)
pairs_sort <- apply(keywords.w[,c(1:2)], 1, str_sort) %>% t
keywords.w <- cbind(pairs_sort, as.double(keywords.w$weight)) %>% as.data.frame()
colnames(keywords.w) <- c("Var1", "Var2", "weight") 
keywords.w$weight <- as.double(keywords.w$weight)

# Sum all values from a same pair
keywords.w <- keywords.w %>% group_by(Var1, Var2) %>% summarise(weight = sum(weight))


# Generate a network to plot
net <- network(keywords.w[, c(1:2)], directed = FALSE)

# Set the edges value
set.edge.attribute(net, "edge", keywords.w$weight)

# Plot the network
set.seed(123)
pdf("Fig1_c.pdf", width=15, height=12)
ggplot(net, aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(color = edge)) +
  geom_nodes(color = "white", size = 0) +
  geom_nodetext(aes(label = vertex.names, size = 5),
                color = "black") +
  scale_color_gradient2(low = "white", high = "purple") +
  guides(size = "none", color = "none") + 
  theme_blank()
dev.off()




