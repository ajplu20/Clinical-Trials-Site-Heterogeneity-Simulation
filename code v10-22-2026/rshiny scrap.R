# Save as app.R and run with: shiny::runApp()

library(shiny)

# --- helpers ---------------------------------------------------------------
parse_num_list <- function(x) {
  # turn comma-separated string into numeric vector; trim spaces
  if (is.null(x) || x == "") return(numeric(0))
  xs <- strsplit(x, ",")[[1]]
  as.numeric(trimws(xs))
}

pretty_valid <- function(ok) if (ok) "✓" else "✗"

# --- UI -------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Power / FP Calculation — Input Prototype"),
  
  sidebarLayout(
    sidebarPanel(
      width = 5,
      h4("Core setup"),
      numericInput("sample_size", "Sample size (per arm or total — your choice):", value = 1000, min = 1, step = 1),
      numericInput("n_sites", "Number of sites:", value = 5, min = 2, step = 1),
      tags$hr(),
      
      h4("Site mortality (control group)"),
      selectInput(
        "mortality_mode", "Input mode:",
        choices = c(
          "Per-site estimates" = "per_site",
          "Overall mean + between-site SD" = "mean_sd"
        ), selected = "per_site"
      ),
      
      # If per-site values provided: comma-separated vector 0..1, length == n_sites
      conditionalPanel(
        condition = "input.mortality_mode == 'per_site'",
        textInput("mortality_list", "Per-site baseline mortalities (comma separated, values in [0,1]):",
                  placeholder = "e.g., 0.10, 0.08, 0.12, 0.09, 0.07")
      ),
      
      # Else ask for overall mean + SD of site-to-site variability
      conditionalPanel(
        condition = "input.mortality_mode == 'mean_sd'",
        numericInput("mortality_mean", "Overall control mortality (0–1):", value = 0.10, min = 0, max = 1, step = 0.001),
        numericInput("mortality_sd", "Between-site SD (>= 0):", value = 0.03, min = 0, step = 0.001)
      ),
      
      tags$hr(),
      h4("Recruitment status"),
      selectInput(
        "recruit_mode", "Input mode:",
        choices = c(
          "Balanced" = "balanced",
          "Expected imbalance (SD around equal shares)" = "imbalance_sd",
          "Custom site proportions" = "custom_props"
        ), selected = "balanced"
      ),
      
      conditionalPanel(
        condition = "input.recruit_mode == 'imbalance_sd'",
        numericInput("recruit_sd", "Recruitment SD (0–1):", value = 0.05, min = 0, max = 1, step = 0.001),
        helpText("SD is relative to equal share 1 / number of sites.")
      ),
      
      conditionalPanel(
        condition = "input.recruit_mode == 'custom_props'",
        textInput("recruit_props", "Per-site recruitment proportions (comma separated, sum to 1):",
                  placeholder = "e.g., 0.30, 0.25, 0.20, 0.15, 0.10")
      ),
      
      tags$hr(),
      h4("Treatment effect across sites"),
      selectInput(
        "te_mode", "Input mode:",
        choices = c(
          "No difference between sites" = "constant",
          "Uniform distribution" = "uniform",
          "Normal distribution" = "normal"
        ), selected = "constant"
      ),
      
      conditionalPanel(
        condition = "input.te_mode == 'constant'",
        numericInput("te_const", "Mean treatment effect:", value = -0.02, step = 0.001)
      ),
      
      conditionalPanel(
        condition = "input.te_mode == 'uniform'",
        numericInput("te_min", "Uniform min:", value = -0.03, min = -1, max = 1, step = 0.001),
        numericInput("te_max", "Uniform max:", value = 0.01,  min = -1, max = 1, step = 0.001)
      ),
      
      conditionalPanel(
        condition = "input.te_mode == 'normal'",
        numericInput("te_mean", "Normal mean (0–1):", value = 0.02, min = 0, max = 1, step = 0.001),
        numericInput("te_sd",   "Normal SD (0–1):",   value = 0.01, min = 0, max = 1, step = 0.001)
      ),
      
      tags$hr(),
      h4("Non-inferiority"),
      checkboxInput("is_ni", "This is a non-inferiority analysis", value = FALSE),
      conditionalPanel(
        condition = "input.is_ni == true",
        numericInput("ni_margin", "Non-inferiority margin (0–1):", value = 0.1, min = 0, max = 1, step = 0.001)
      )
    ),
    
    mainPanel(
      width = 7,
      h4("Validation status"),
      verbatimTextOutput("status"),
      uiOutput("calcUI"),
      conditionalPanel(
        condition = "input.calc > 0",
        h4("Simulation progress"),
        uiOutput("progressBar"),
        br(), br(),
        h4("Results (placeholders)"),
        tags$p(strong("Power:"), span("90%")),
        tags$p(strong("False positive rate:"), span("5%")),
        tags$hr()
      )
    )
  )
)

# --- Server ----------------------------------------------------------------
server <- function(input, output, session) {
  # Collect parameters as a reactive list (keeps dropdown choices as factors if desired)
  params <- reactive({
    list(
      sample_size = input$sample_size,
      n_sites     = input$n_sites,
      mortality_mode = factor(input$mortality_mode, levels = c("per_site", "mean_sd")),
      mortality_list = parse_num_list(input$mortality_list),
      mortality_mean = input$mortality_mean,
      mortality_sd   = input$mortality_sd,
      recruit_mode   = factor(input$recruit_mode, levels = c("balanced", "imbalance_sd", "custom_props")),
      recruit_sd     = input$recruit_sd,
      recruit_props  = parse_num_list(input$recruit_props),
      te_mode        = factor(input$te_mode, levels = c("constant", "uniform", "normal")),
      te_const       = input$te_const,
      te_min         = input$te_min,
      te_max         = input$te_max,
      te_mean        = input$te_mean,
      te_sd          = input$te_sd,
      is_ni          = isTRUE(input$is_ni),
      ni_margin      = input$ni_margin
    )
  })

# Validation logic (returns character vector of messages)
validate_inputs <- reactive({
  p <- params()
  msgs <- c()
  
  # basics
  if (is.na(p$sample_size) || p$sample_size < 1) msgs <- c(msgs, "Sample size must be >= 1.")
  if (is.na(p$n_sites) || p$n_sites < 2) msgs <- c(msgs, "Number of sites must be >= 2.")
  
  # mortality
  if (p$mortality_mode == "per_site") {
    v <- p$mortality_list
    if (length(v) != p$n_sites) msgs <- c(msgs, sprintf("Mortality list length (%d) must equal number of sites (%d).", length(v), p$n_sites))
    if (any(is.na(v))) msgs <- c(msgs, "Mortality list contains non-numeric values.")
    if (length(v) > 0 && (min(v) < 0 || max(v) > 1)) msgs <- c(msgs, "Mortality values must be in [0,1].")
  } else {
    if (is.na(p$mortality_mean) || p$mortality_mean < 0 || p$mortality_mean > 1) msgs <- c(msgs, "Overall mortality must be in [0,1].")
    if (is.na(p$mortality_sd) || p$mortality_sd < 0) msgs <- c(msgs, "Between-site SD must be >= 0.")
  }
  
  # recruitment
  if (p$recruit_mode == "imbalance_sd") {
    if (is.na(p$recruit_sd) || p$recruit_sd < 0 || p$recruit_sd > 1) msgs <- c(msgs, "Recruitment SD must be in [0,1].")
  } else if (p$recruit_mode == "custom_props") {
    r <- p$recruit_props
    if (length(r) != p$n_sites) msgs <- c(msgs, sprintf("Recruitment vector length (%d) must equal number of sites (%d).", length(r), p$n_sites))
    if (any(is.na(r))) msgs <- c(msgs, "Recruitment vector contains non-numeric values.")
    if (length(r) > 0 && (min(r) < 0 || max(r) > 1)) msgs <- c(msgs, "Recruitment proportions must be in [0,1].")
    if (length(r) > 0 && abs(sum(r) - 1) > 1e-6) msgs <- c(msgs, sprintf("Recruitment proportions must sum to 1 (currently %.6f).", sum(r)))
  }
  
  # treatment effect
  if (p$te_mode == "uniform") {
    if (!is.na(p$te_min) && !is.na(p$te_max) && p$te_min > p$te_max) msgs <- c(msgs, "Uniform min must be <= max.")
    if (is.na(p$te_min) || p$te_min < -1 || p$te_min > 1) msgs <- c(msgs, "Uniform min must be in [-1,1].")
    if (is.na(p$te_max) || p$te_max < -1 || p$te_max > 1) msgs <- c(msgs, "Uniform max must be in [-1,1].")
  } else if (p$te_mode == "normal") {
    if (is.na(p$te_mean) || p$te_mean < 0 || p$te_mean > 1) msgs <- c(msgs, "Normal mean must be in [0,1].")
    if (is.na(p$te_sd)   || p$te_sd   < 0 || p$te_sd   > 1) msgs <- c(msgs, "Normal SD must be in [0,1].")
  }
  
  # NI
  if (p$is_ni) {
    if (is.na(p$ni_margin) || p$ni_margin < 0 || p$ni_margin > 1) msgs <- c(msgs, "NI margin must be in [0,1].")
  }
  
  msgs
})

# Show validation text
output$status <- renderText({
  msgs <- validate_inputs()
  if (length(msgs) == 0) {
    "All inputs look valid. Click 'Calculate' to proceed."
  } else {
    paste("Please address the following:", paste(paste0("- ", msgs), collapse = "
"), sep = "
")
  }
})

# Calculate button: rendered only when there are no validation errors
output$calcUI <- renderUI({
  if (length(validate_inputs()) == 0) {
    actionButton("calc", "Calculate", class = "btn btn-primary")
  } else {
    tags$div(class = "text-muted", "Fix inputs to enable the Calculate button.")
  }
})

# Static 50% progress bar (only appears after Calculate is pressed)
output$progressBar <- renderUI({
  tags$div(
    style = "width:100%; height:24px; background:#eee; border-radius:12px; overflow:hidden; box-shadow: inset 0 0 3px rgba(0,0,0,0.2);",
    tags$div(style = "height:100%; width:50%; background:#7fc97f; transition: width 0.3s;")
  )
})

}

shinyApp(ui, server)
