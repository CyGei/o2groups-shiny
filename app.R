library(shiny)
library(shinythemes)
library(shinycssloaders)
library(plotly)
library(tidyverse)
library(magrittr)
library(here)
library(o2groups)
library(rhandsontable)

linebreaks <- function(n) {
  HTML(strrep(br(), n))
}

setwd(here::here())
source("helpers.R")

ui <- fluidPage(
  theme = shinytheme("slate"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "n_groups",
        "Number of groups",
        min = 1,
        max = 10,
        value = 4
      ),
      linebreaks(2),
      numericInput("mean_gen_time", "Mean Generation Time", value = 6.2),
      numericInput("sd_gen_time", "Standard Deviation Generation Time", value = 1.5),
      numericInput("mean_incubation", "Mean Incubation Period", value = 5.1),
      numericInput("sd_incubation", "Standard Deviation Incubation Period", value = 1.5),
      linebreaks(2),
      rHandsontableOutput("param_table"),
      actionButton("run", "Run simulation"),
      downloadButton("downloadData", "Download Data"),
      downloadButton("downloadStats", "Download Stats"),
      linebreaks(2)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Tree Plot", plotlyOutput("p_tree", height = "100vh") %>% shinycssloaders::withSpinner(color = "#0dc5c1")),
        tabPanel("Incidence Plot", plotOutput("p_incidence", height = "100vh")),
        tabPanel("Matrices Plot", plotOutput("p_matrices", height = "100vh")),
        tabPanel("Results Plot", plotOutput("p_results", height = "100vh"))
      )
    )
  )
)


server <- function(input, output) {
  # Create a reactive dataset for the parameter table
  param_table_data <- reactive({
    n_groups <- input$n_groups
    data <- data.frame(
      Name = LETTERS[1:n_groups],
      Size = rep("", n_groups),
      Delta = rep("", n_groups),
      R0 = rep("", n_groups),
      Introductions = rep("", n_groups)
    )
    data
  })
  
  # Render the parameter table
  output$param_table <- renderRHandsontable({
    rhandsontable(param_table_data())
  })
  
  # Generate data
  out <- eventReactive(input$run, {
    param_data <- hot_to_r(input$param_table)
    n_groups <- nrow(param_data)
    name <- param_data$Name
    size <- as.numeric(param_data$Size)
    delta <- as.numeric(param_data$Delta)
    r0 <- as.numeric(param_data$R0)
    intro_n <- as.numeric(param_data$Introductions)
    
    generation_time <- c(input$mean_gen_time, input$sd_gen_time)
    incubation_period <- c(input$mean_incubation, input$sd_incubation)
    
    # run simulation
    set.seed(123)
    out <- o2groups::simulate_groups_furrr(
      n_simulations = 100,
      duration = 100,
      n_groups = n_groups,
      name = name,
      size = size,
      delta = delta,
      r0 = r0,
      intro_n = intro_n,
      generation_time = simulacr::make_disc_gamma(generation_time[1], generation_time[2])$d(0:100),
      incubation_period = simulacr::make_disc_gamma(incubation_period[1], incubation_period[2])$r(1000)
    )
    
    param_info <- list(
      n_groups = n_groups,
      name = name,
      size = size,
      delta = delta,
      r0 = r0,
      intro_n = intro_n
    )
    
    
    return(list(
      data = out[[1]],
      stats = out[[2]],
      param_info = param_info
    ))
  })
  
  
  # Plot
  output$p_tree <- renderPlotly({
    data <- out()$data
    p_tree <- plot_tree(subset(data, simulation == 1),
                        pal = scales::hue_pal()(out()$param_info$n_groups)
    )
    p_tree
  })
  
  output$p_incidence <- renderPlot({
    data <- out()$data
    p_incidence <- plot_stats(data)
    p_incidence
  })
  
  output$p_matrices <- renderPlot({
    data <- out()$data
    param_info <- out()$param_info
    p_Ri <- p_Ri(data, param_info)
    p_mixing <- p_mixing(data, param_info)
    p_matrices <- patchwork::wrap_plots(p_Ri, p_mixing, nrow = 2)
    p_matrices
  })
  
  output$p_results <- renderPlot({
    data <- out()$data
    param_info <- out()$param_info
    p_results <- p_delta(data, param_info)
    p_results
  })
  
  
  # Download data:
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(out()$data, file, row.names = FALSE)
    }
  )
  
  # Download Stats:
  output$downloadStats <- downloadHandler(
    filename = function() {
      paste("stats-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(out()$stats, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
