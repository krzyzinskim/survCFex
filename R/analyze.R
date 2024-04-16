#' @export
analyze <- function(counterfactual_explanations) {
  # Check data structure
  if (!is(counterfactual_explanations, "counterfactual_explanations") |
      !is.data.frame(counterfactual_explanations$original_observation) |
      !is.data.frame(counterfactual_explanations$counterfactual_examples) |
      !is.data.frame(counterfactual_explanations$objective_values)) {
    stop("Invalid data structure. Please provide an object of class 'counterfactual_explanations' with the following components: 'original_observation', 'counterfactual_examples', 'objective_values'.")
  }

  loss_functions <- names(counterfactual_explanations$objective_values)
  type <- class(counterfactual_explanations)[2]

  # Define UI elements
  ui <- fluidPage(
    titlePanel("Counterfactual Analysis"),
    sidebarLayout(
      sidebarPanel(
        h4("Filter counterfactuals by objective values"),
        h6("Choose the maximum value for each loss function."),
        p("Note that all objective values are minimized."),
        lapply(loss_functions, function(loss) {
          sliderInput(
            paste0(loss, "_slider"),
            h5(loss),
            min = 0,
            max = ceiling(max(counterfactual_explanations$objective_values[, loss]) * 100) / 100,
            value = ceiling(max(counterfactual_explanations$objective_values[, loss]) * 100) / 100,
            step = ifelse(loss == "sparsity", 1, 0.01),
            round = ifelse(loss == "sparsity", TRUE, 2)
          )
        }),
        width = 3
      ),
      mainPanel(
        # original observation
        h3("Original Observation"),
        DT::dataTableOutput("original_observation"),
        # counterfactuals
        h3("Counterfactuals"),
        h5("There are",
           textOutput("counterfactuals_count", inline = TRUE),
           "counterfactuals fulfilling the criteria."),
        tabsetPanel(
          tabPanel("Table",
                   DT::dataTableOutput("counterfactuals")),
          tabPanel("Visualizations",
                   radioButtons(
                     "function_type",
                     label = "Choose prediction function:",
                     choices = c("Survival Function" = "survival",
                                 "Cumulative Hazard Function" = "chf"),
                     selected = ifelse(type == "multiobjective_counterfactuals", "survival", "chf")
                   ),
                   plotly::plotlyOutput("counterfactuals_plot", height="400px"),
                   plotOutput("parallel_coordinates", height="400px"),
                   plotOutput("changes_frequency", height="300px"))
        )
      )
    )
  )

  # Define server logic
  server <- function(input, output) {
    filtered_data <- reactive({
      filtered_population <- cbind(counterfactual_explanations$counterfactual_examples, counterfactual_explanations$objective_values)
      for (loss in loss_functions) {
        filtered_population <- filtered_population[filtered_population[, loss] <=
                                                     input[[paste0(loss, "_slider")]], ]
      }
      return(filtered_population)
    })

    filtered_predictions <- reactive({
      which_rows <- rownames(counterfactual_explanations$counterfactual_examples) %in% rownames(filtered_data())
      filtered_preds <- counterfactual_explanations$predictions[which_rows,]
      return(filtered_preds)
    })

    output$original_observation <- DT::renderDataTable({
      DT::datatable(
        counterfactual_explanations$original_observation,
        options = list(ordering = FALSE, paging = FALSE, dom = 't'),
        rownames = FALSE
      ) |>
        DT::formatStyle(
          columns = 1:ncol(counterfactual_explanations$original_observation),
          backgroundColor = 'ghostwhite'
        )
    })

    output$counterfactuals_plot <- renderPlotly({
      plot_counterfactual_predictions(counterfactual_explanations,
                                      filtered_predictions(),
                                      filtered_data(),
                                      input$function_type)
    })

    output$parallel_coordinates <- renderPlot({
      plot_parallel_coordinates(counterfactual_explanations, filtered_data())
    })

    output$changes_frequency <- renderPlot({
      plot_changes_frequency(counterfactual_explanations, filtered_data())
    })

    output$counterfactuals <- DT::renderDataTable({
      DT::datatable(
        filtered_data(),
        options = list(ordering = TRUE, paging = FALSE, dom = 't'),
        rownames = TRUE
      ) |>
        DT::formatRound(
          columns = c("validity", "similarity", "plausibility"),
          digits = 3
        ) |>
        DT::formatRound(
          columns = "sparsity",
          digits = 0
        ) |>
        DT::formatStyle(
          columns = 1:ncol(counterfactual_explanations$counterfactual_examples),
          backgroundColor = 'ghostwhite'
        ) |>
        DT::formatStyle(
          columns = (ncol(counterfactual_explanations$counterfactual_examples) + 1):ncol(filtered_data()),
          backgroundColor = 'lavenderblush',
        )
    })

    output$counterfactuals_count <- renderText({
      nrow(filtered_data())
    })
  }

  # Run the app
  shinyApp(ui = ui, server = server)
}
