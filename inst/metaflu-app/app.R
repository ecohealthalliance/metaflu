#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(metaflu)
library(ggplot2)
library(dplyr)
library(doMC)
library(tidyr)
library(purrr)
library(gridExtra)
library(abind)
library(doRNG)
library(igraph)
library(threejs)
library(shiny)
library(metaflu)
library(ggplot2)
library(dplyr)
library(purrr)
library(gridExtra)
library(abind)
library(grid)
set.seed(123)

#Set number of farms and ~ number of chickens per farm
farm_number <- 100
farm_size <- 40



metaflu_sim <- function(sims, cull_time){
  parms = list(
    beta = 1.44456,   #contact rate for direct transmission
    gamma = 0.167,  #recovery rate
    mu = 0,         #base mortality rate
    alpha = 0.4,      #disease mortality rate
    phi = 0,  #infectiousness of environmental virions
    eta = 0,     #degradation rate of environmental virions
    nu =  0.00,    #uptake rate of environmental virion
    sigma = 0,      #virion shedding rate
    omega = 0.03,   #movement rate
    rho = 0.85256,        #contact  nonlinearity 0=dens-dependent, 1=freq-dependent
    lambda = 0,     #force of infection from external sources
    tau_crit = 5,   #critical suveillance time
    I_crit = 1,     #threshold for reporting
    pi_report = 0.1, #reporting probability
    pi_detect = 0.9, #detection probability
    cull_time = cull_time,   #time to detect
    network_type = "smallworld",
    network_parms = list(dim = 1, size = farm_number, nei = 2.33, p = 0.0596, multiple = FALSE, loops = FALSE),
    stochastic_network = TRUE
  )

  g_list <- parallel::mclapply(seq_len(sims), function(y){
    patches <- grow_patches_clustered(basic_patches(40,100))
    i_patches <- seed_initial_infection(patches)
    result <- mf_sim(init = i_patches, parameters = parms, times=1:365, n_sims = 1)
    return(result)
  }, mc.cores = parallel::detectCores())
  g_results <- do.call("abind", g_list)

  return(list(g_results, g_list))

}


#################### SHINY ########################

ui <- fluidPage(
  titlePanel("Metaflu Clustered Growth Simulator"),
  h5("This application allows for the side-by-side comparisons of simulations for different values of the culling time parameter."),
  HTML("<br>"),
  fluidRow(column(6,
                  #sims
                  sliderInput(inputId = "sims",
                              label = "Number of Simulations",
                              value = 10, min = 10, max = 500, step = 10),

                  #cull-time
                  sliderInput(inputId = "cull_time",
                              label = "Select Time to Culling in Days",
                              value = 1, min = 1, max = 20),

                  actionButton(inputId = "update",
                               label = "Simulate"),

                  plotOutput("graph_panel", width = "100%", height = "500px"),
                  actionButton(inputId = "reset_3d",
                               label = "Animate Random Run"),
                  scatterplotThreeOutput("threed", width = "100%", height = "500px")


  ),

  column(6,
         #sims
         sliderInput(inputId = "sims_2",
                     label = "Number of Simulations",
                     value = 10, min = 10, max = 500, step = 10),

         #cull-time
         sliderInput(inputId = "cull_time_2",
                     label = "Select Time to Culling in Days",
                     value = 1, min = 1, max = 20),

         actionButton(inputId = "update_2",
                      label = "Simulate"),
         plotOutput("graph_panel_2", width = "100%", height = "500px"),
         actionButton(inputId = "reset_3d_2",
                      label = "Animate Random Run"),
         scatterplotThreeOutput("threed_2", width = "100%", height = "500px")

  )
  )


)
server <- function(input, output){
  #data <- eventReactive(input$update, {
  # create_graphjs(g_list[[isolate({input$num})]], isolate({input$num}))
  #})


  data <- eventReactive(input$update, {
    metaflu_sim(input$sims, input$cull_time)
  })

  data2 <- eventReactive(input$update_2, {
    metaflu_sim(input$sims_2, input$cull_time_2)
  })

  rando_graph <- eventReactive(input$reset_3d, {
    sample.int(length(data()[[2]]),1)
  })

  rando_graph_2 <- eventReactive(input$reset_3d_2, {
    sample.int(length(data2()[[2]]),1)
  })

  output$graph_panel <-
    renderPlot({
      create_graph_panel(data()[[1]], paste(input$sims, "Simulations with Cull Time of", input$cull_time, "Days"))
    })

  output$threed <- renderScatterplotThree({
    if (is.null(rando_graph)) return()
    create_graphjs(data()[[2]][[rando_graph()]], rando_graph())
  })

  output$graph_panel_2 <-

    renderPlot({
      create_graph_panel(data2()[[1]], paste(input$sims_2, "Simulations with Cull Time of", input$cull_time_2, "Days"))
    })

  output$threed_2 <- renderScatterplotThree({
    if (is.null(rando_graph_2)) return()
    create_graphjs(data()[[2]][[rando_graph_2()]], rando_graph_2())
  })

}



shinyApp(ui = ui, server = server)
