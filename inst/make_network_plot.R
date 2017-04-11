library(igraph)
library(ggplot2)
library(ggraph)
library(ggforce)
library(magrittr)
#library(tidyverse)
#net_adj <- readRDS("epinetwork.rds")$network
#net <- graph_from_adjacency_matrix(net_adj > 0) # Use this to generate a graph
                                                 # from simulation network outputs

set.seed(20)
net <- sample_smallworld(dim = 1, size = 200, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE)
net <- set_vertex_attr(net, "big", value=as.factor(sample(c(rep(0, 180), rep(1, 20)), 200)))
net <- set_vertex_attr(net, "big2", value=as.factor(c(rep(0, 40), rep(1, 20), rep(0, 140))))
net <- set_vertex_attr(net, "label", value=as.character(1:200))

mylay <- create_layout(net, layout = "igraph", algorithm="kk")
base <- ggraph(net,'manual', node.positions= mylay[,c("x", "y")]) +
  geom_edge_arc0(curvature=0.1, edge_width=0.5) +
  geom_node_point(size=5, pch=21, colour="black", fill="red", stroke=1.5) +
  theme_void()
ggsave(base, filename = "base_network.png", height=5, width=5)

growth <- ggraph(net,'manual', node.positions= mylay[,c("x", "y")]) +
  geom_edge_arc0(curvature=0.1, edge_width=0.5) +
  geom_node_point(mapping=aes(size=big, fill=as.factor(big)), pch=21, colour="black", stroke=1.5) +
  scale_size_continuous(range = c(5,15)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_void() +
  theme(legend.position="none")
ggsave(growth, filename = "growth_network.png", height=5, width=5)

growth_local <- ggraph(net,'manual', node.positions= mylay[,c("x", "y")]) +
  geom_edge_arc0(curvature=0.1, edge_width=0.5) +
  geom_node_point(mapping=aes(size=big2, fill=as.factor(big2)), pch=21, colour="black", stroke=1.5) +
  scale_size_continuous(range = c(5,15)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_void() +
  theme(legend.position="none")
ggsave(growth_local, filename = "growth_local_network.png", height=5, width=5)

net2 <- subgraph.edges(net, 1:266, delete.vertices = FALSE)
pruned <- ggraph(net2,'manual', node.positions= mylay[,c("x", "y")]) +
  geom_edge_arc0(curvature=0.1, edge_width=0.5) +
  geom_node_point(size=5, pch=21, colour="black", fill="red", stroke=1.5) +
  theme_void()
ggsave(pruned, filename = "pruned_network.png", height=5, width=5)

