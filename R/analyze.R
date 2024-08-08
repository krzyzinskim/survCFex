#' @title Analyze counterfactual explanations using a Shiny app
#'
#' @param counterfactual_explanations An object of class `counterfactual_explanations` containing the original observation, generated counterfactuals, and objective values.
#' @param digits Number of digits to round the numerical values to.
#'
#' @import shiny
#' @export

analyze <- function(counterfactual_explanations, digits = 2) {
  # Check data structure
  if (!inherits(counterfactual_explanations,
                "counterfactual_explanations") |
      !is.data.frame(counterfactual_explanations$original_observation) |
      !is.data.frame(counterfactual_explanations$counterfactual_examples) |
      !is.data.frame(counterfactual_explanations$objective_values)) {
    stop(
      "Invalid data structure. Please provide an object of class 'counterfactual_explanations' with the following components: 'original_observation', 'counterfactual_examples', 'objective_values'."
    )
  }
  loss_functions <-
    names(counterfactual_explanations$objective_values)
  type <- class(counterfactual_explanations)[2]

  num_ids <-
    which(sapply(counterfactual_explanations$counterfactual_examples, function(x) {
      is.numeric(x)
    }))
  to_round_ids <-
    num_ids[sapply(num_ids, function(i) {
      any(counterfactual_explanations$counterfactual_examples[, i] %% 1 != 0)
    })]

  # Define UI elements
  ui <- fluidPage(titlePanel("Counterfactual Analysis"),
                  sidebarLayout(
                    sidebarPanel(
                      h3("Filter counterfactuals by objective values"),
                      h5("Choose the maximum value for each loss function."),
                      h6("Note that all objective values are minimized."),
                      lapply(loss_functions, function(loss) {
                        sliderInput(
                          paste0(loss, "_slider"),
                          h5(paste(loss, "loss")),
                          min = 0,
                          max = ceiling(
                            max(counterfactual_explanations$objective_values[, loss]) * 10 ^ digits
                          ) / 10 ^ digits,
                          value = ceiling(
                            max(counterfactual_explanations$objective_values[, loss]) * 10 ^ digits
                          ) / 10 ^ digits,
                          step = ifelse(loss == "sparsity", 1, 10 ^ (-digits)),
                          round = ifelse(loss == "sparsity", TRUE, digits)
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
                      h5(
                        "There are",
                        textOutput("counterfactuals_count", inline = TRUE),
                        "counterfactuals fulfilling the criteria."
                      ),
                      tabsetPanel(
                        tabPanel("Table",
                                 DT::dataTableOutput("counterfactuals")),
                        tabPanel(
                          "Visualizations",
                          radioButtons(
                            "function_type",
                            label = "Choose prediction function:",
                            choices = c(
                              "Survival Function" = "survival",
                              "Cumulative Hazard Function" = "chf"
                            ),
                            selected = ifelse(type == "multiobjective_counterfactuals", "survival", "chf")
                          ),
                          plotlyOutput("counterfactuals_plot", height = "400px"),
                          plotOutput("parallel_coordinates", height = "400px"),
                          plotOutput("changes_frequency", height = "300px")
                        )
                      )
                    )
                  ))

  # Define server logic
  server <- function(input, output) {
    filtered_data <- reactive({
      filtered_population <-
        cbind(
          counterfactual_explanations$counterfactual_examples,
          counterfactual_explanations$objective_values
        )
      for (loss in loss_functions) {
        filtered_population <-
          filtered_population[filtered_population[, loss] <=
                                input[[paste0(loss, "_slider")]], , drop =
                                FALSE]
      }
      return(filtered_population)
    })

    filtered_predictions <- reactive({
      which_rows <-
        rownames(counterfactual_explanations$counterfactual_examples) %in% rownames(filtered_data())
      filtered_preds <-
        counterfactual_explanations$predictions[which_rows, , drop = FALSE]
      return(filtered_preds)
    })

    output$original_observation <- DT::renderDataTable({
      DT::datatable(
        counterfactual_explanations$original_observation,
        options = list(
          ordering = FALSE,
          paging = FALSE,
          dom = 't'
        ),
        rownames = FALSE
      ) |>
        DT::formatRound(columns = to_round_ids,
                        digits = digits) |>
        DT::formatStyle(
          columns = 1:ncol(counterfactual_explanations$original_observation),
          backgroundColor = 'ghostwhite',
        )
    })

    output$counterfactuals_plot <- renderPlotly({
      plot_counterfactual_predictions(
        counterfactual_explanations,
        filtered_predictions(),
        filtered_data(),
        input$function_type
      )
    })

    output$parallel_coordinates <- renderPlot({
      plot_parallel_coordinates(counterfactual_explanations, filtered_data())
    },
    res = 90)

    output$changes_frequency <- renderPlot({
      plot_changes_frequency(counterfactual_explanations, filtered_data())
    },
    res = 90)

    output$counterfactuals <- DT::renderDataTable({
      DT::datatable(
        filtered_data(),
        options = list(
          ordering = TRUE,
          paging = FALSE,
          dom = 't'
        ),
        rownames = TRUE
      ) |>
        DT::formatRound(
          columns = c("validity", "similarity", "plausibility"),
          digits = digits
        ) |>
        DT::formatRound(columns = "sparsity",
                        digits = 0) |>
        DT::formatRound(columns = to_round_ids,
                        digits = digits) |>
        DT::formatStyle(
          columns = 1:ncol(counterfactual_explanations$counterfactual_examples),
          backgroundColor = 'ghostwhite',
        ) |>
        DT::formatStyle(
          columns = (
            ncol(counterfactual_explanations$counterfactual_examples) + 1
          ):ncol(filtered_data()),
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
