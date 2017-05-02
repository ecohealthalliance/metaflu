#' @export
#' @importFrom igraph graph_from_adjacency_matrix layout_nicely
#' @importFrom threejs graphjs
#'
create_graphjs <- function(sim, num){
  net_adj <- attr(sim, "network")$`1`
  net <- graph_from_adjacency_matrix(net_adj > 0)
  graph_layout <- layout_nicely(net,3)
  sim_dur <- get_duration_array(sim)$duration
  graph <- lapply(seq_len(sim_dur + 1), function(x){
    S <- sim["S",,x,1]
    I <- sim["I",,x,1]
    C <- sim["C",,x,1]
    net <- graph_from_adjacency_matrix(net_adj > 0)
    V(net)$size <- log2(S)
    V(net)$color <- "grey"
    inf <- which(I > 0)
    cul <- which(C > 1)
    V(net)[inf]$color <- "red"
    V(net)[cul]$color <- "blue"
    return(net)
  })
  graphjs(graph, main = paste("Simulation:", num, "Day ",1:365), fpl = 60, layout = graph_layout)
}

#' @export
#' @importFrom plyr adply
#' @importFrom dplyr %>% rename group_by mutate filter select
#' @importFrom tidyr spread
#' @importFrom igraph graph_from_adjacency_matrix layout_with_graphopt
#' @importFrom ggraph ggraph geom_edge_arc geom_node_point
#' @importFrom gganimate gganimate
#' @importFrom animation ani.options

create_gif <- function(list_el, fname, g_status = FALSE, title = ""){
  sim_dur <- get_duration_array(list_el)$duration
  simulation <- list_el[,,1:(sim_dur+1),1]
  sim_net <- attr(list_el, "network")$`1`
  sim_df <- plyr::adply(simulation, c(1,2,3)) %>%
    dplyr::rename(compartment = X1, patch = X2, time = X3, population = V1) %>%
    tidyr::spread(compartment, population)
  raw_net <- graph_from_adjacency_matrix(sim_net > 0, "undirected")

  layout <- layout_nicely(raw_net,2) %>% as.data.frame() %>% rename(x=V1, y=V2)

  sim1 <- sim_df %>%
    group_by(time) %>%
    mutate(x = as.numeric(layout$x), y = as.numeric(layout$y)) %>%
    group_by()

  sim1 <- sim1 %>%
    group_by(patch,time) %>%
    mutate(status = if_else(I > 0, "infected", if_else(C > 1, "culled", "normal")))

  original_s <- sim1 %>%
    filter(time == 1) %>%
    select(patch, original_S = S, blargo = time)

  sim1 <- left_join(sim1, original_s, by = "patch") %>%
    select(-blargo)

  sim_end = filter(sim1, time==1) # A single time point for use for background

  cols <- c("normal" = "black", "infected" = "red", "culled" = "blue")

  if(g_status == TRUE){
    size_range <- c(4,12)
  }else{
    size_range <- c(5,7)
  }

  plot <- ggraph(raw_net, 'manual', node.positions=layout) +
    geom_edge_arc(curvature = 0.1, edge_colour='grey40') +
    geom_node_point(data=sim_end, mapping=aes(x=x, y=y,size=original_S), #white bg for the nodes
                    shape=21, stroke=1,          #to hide the edges beind
                    color="white", fill="white") +       #them
    geom_node_point(data=sim1,
                    mapping=aes(x=x, y=y, fill= status, alpha=(S/original_S)*100, frame=time, size = original_S),
                    shape=21, stroke=1, color="grey10") +
    scale_fill_manual(values = cols) +
    scale_alpha_continuous(range=c(0.5,1)) +
    scale_size_continuous(range=size_range) +
    labs(title = title) +
    ggforce::theme_no_axes() +
    theme(legend.position="none",
          panel.border = element_blank())

  animation::ani.options(interval=0.25)
  gganimate(plot, paste0(fname,'.gif'), title_frame = FALSE)
}


#' @export
#' @importFrom gridExtra grid.arrange

create_graph_panel <- function(result_array, title = ""){
  gdurations <- get_duration_array(result_array)
  gexposure <- get_exposure(result_array)
  gs.df <- get_all_sims("S", result_array)
  gi.df <- get_all_sims("I", result_array)
  gr.df <- get_all_sims("R", result_array)
  gfarm_no <- get_number_farms(result_array)
  gloss <- get_proportion_loss(result_array)

  sir.graph <- ggplot() +
    geom_line(data = gs.df, aes(x = time, y = pop, group = sim), alpha = 0.05, color = "blue") +
    geom_line(data = gi.df, aes(x = time, y = pop, group = sim), alpha = 0.05, color = "red") +
    geom_line(data = gr.df, aes(x = time, y = pop, group = sim), alpha = 0.05, color = "green") +
    labs(title = paste("S I R"), x = "Time", y = "Population") +
    scale_colour_manual(name = "Compartment", values=c(S = "blue", I = "red", R = "green")) +
    theme_minimal()

  exposure.graph <- ggplot() +
    geom_histogram(data = gexposure, aes(x = inf_exp), bins = 100) +
    labs(title = "Exposure", x = "Infected Chicken-Days", y = "Number of Simulations") +
    theme_minimal()

  farm.graph <- ggplot() +
    geom_histogram(data = gfarm_no, aes(x = num_farms), bins = 100) +
    labs(title = "Farm Spread", x = "Farms Infected", y = "Number of Simulations") +
    coord_flip() +
    theme_minimal()

  duration.graph <- ggplot() +
    geom_histogram(data = gdurations, aes(x = duration), bins = 100) +
    theme_minimal() +
    labs(title = "Duration of Epidemic", x = "Duration in Days", y = "Number of Simulations")

  lay <- rbind(c(1,1,1,4),
               c(1,1,1,4),
               c(2,2,3,3))

  grid.arrange(sir.graph,exposure.graph,duration.graph,farm.graph, layout_matrix = lay, top=textGrob(title, gp = gpar(fontsize = 16)))

}


