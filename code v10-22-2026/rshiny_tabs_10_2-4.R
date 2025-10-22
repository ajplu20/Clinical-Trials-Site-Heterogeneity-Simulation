# Save as app.R and run with: shiny::runApp()

library(shiny)
library(ggplot2)

# ---------- helpers ----------
`%||%` <- function(a, b) if (!is.null(a)) a else b
nz <- function(x, na = NA_real_) if (is.null(x)) na else x

parse_num_list <- function(x) {
  if (is.null(x) || x == "") return(numeric(0))
  xs <- strsplit(x, ",")[[1]]
  as.numeric(trimws(xs))
}

parse_int_list <- function(x) {
  vals <- parse_num_list(x)
  vals <- vals[!is.na(vals)]
  as.integer(vals)
}

# ---------- Reusable UI blocks (prefix-aware) ----------
site_mortality_block <- function(prefix = "") {
  tagList(
    h4("Site baseline mortality (control group)"),
    selectInput(
      paste0(prefix, "mortality_mode"), "How to set baseline mortality:",
      choices = c("Per-site estimates" = "per_site",
                  "Uniform across a range" = "uniform_range",
                  "Overall mean + between-site SD" = "mean_sd"),
      selected = "per_site"
    ),
    
    # mean/sd mode
    conditionalPanel(
      condition = sprintf("input['%smortality_mode'] == 'mean_sd'", prefix),
      sliderInput(paste0(prefix, "mortality_mean"),
                  "Average control mortality:",
                  min = 0, max = 1, value = 0.10, step = 0.005),
      helpText("Sets the overall average mortality across all sites."),
      sliderInput(paste0(prefix, "mortality_sd"),
                  "Between-site variation (SD):",
                  min = 0, max = 1, value = 0.03, step = 0.005),
      helpText("Controls how different sites are from each other around the average.")
    ),
    
    # per-site mode: dynamic rows for each site count
    conditionalPanel(
      condition = sprintf("input['%smortality_mode'] == 'per_site'", prefix),
      uiOutput(paste0(prefix, "mortality_lists_container")),
      helpText("Enter the baseline mortality for each site directly.")
    ),
    
    # uniform range mode: range slider min/max
    conditionalPanel(
      condition = sprintf("input['%smortality_mode'] == 'uniform_range'", prefix),
      sliderInput(paste0(prefix, "mortality_range"),
                  "Mortality range (min to max):",
                  min = 0, max = 1, value = c(0.05, 0.15), step = 0.005),
      helpText("Sites will be spread evenly between the minimum and maximum mortality values.")
    )
  )
}

recruitment_block <- function(prefix = "") {
  tagList(
    h4("Recruitment status"),
    selectInput(
      paste0(prefix, "recruit_mode"), "How to set recruitment across sites:",
      choices = c("Balanced (equal at all sites)" = "balanced",
                  "Expected imbalance (around equal shares)" = "imbalance_sd",
                  "Custom proportions" = "custom_props"),
      selected = "balanced"
    ),
    
    # imbalance_sd mode
    conditionalPanel(
      condition = sprintf("input['%srecruit_mode'] == 'imbalance_sd'", prefix),
      sliderInput(paste0(prefix, "recruit_sd"),
                  "Recruitment variation (SD):",
                  min = 0, max = 1, value = 0.05, step = 0.005),
      helpText("Larger means recruitment is more uneven across sites compared to equal share (1/[site count] per site).")
    ),
    
    # custom proportions mode
    conditionalPanel(
      condition = sprintf("input['%srecruit_mode'] == 'custom_props'", prefix),
      uiOutput(paste0(prefix, "recruit_props_container")),
      helpText("Enter the exact recruitment proportion for each site (must sum to 1).")
    )
  )
}

treatment_effect_block <- function(prefix = "") {
  tagList(
    h4("Treatment effect across sites"),
    helpText("How much lower the treatment group’s mortality is compared to the control group"),
    selectInput(
      paste0(prefix, "te_mode"), "Input mode:",
      choices = c("Constant (same at all sites)" = "constant",
                  "Uniform (vary randomly from a range)" = "uniform",
                  "Normal (vary randomly around a mean)" = "normal"),
      selected = "constant"
    ),
    
    # constant effect mode
    conditionalPanel(
      condition = sprintf("input['%ste_mode'] == 'constant'", prefix),
      sliderInput(paste0(prefix, "te_const"),
                  "Average treatment effect:",
                  min = -1, max = 1, value = 0.10, step = 0.01),
      helpText("Sets one fixed effect size used for every site.")
    ),
    
    # uniform effect mode (distinct IDs)
    conditionalPanel(
      condition = sprintf("input['%ste_mode'] == 'uniform'", prefix),
      sliderInput(paste0(prefix, "te_center_uniform"),
                  "Center of the range:",
                  min = -1, max = 1, value = 0.10, step = 0.01),
      helpText("The midpoint of the treatment effect range across sites."),
      sliderInput(paste0(prefix, "te_width_uniform"),
                  "Width of the range:",
                  min = 0, max = 1, value = 0.20, step = 0.01),
      helpText("Controls how wide the spread of the range is around the center value.")
    ),
    
    # normal effect mode (distinct IDs)
    conditionalPanel(
      condition = sprintf("input['%ste_mode'] == 'normal'", prefix),
      sliderInput(paste0(prefix, "te_mean_normal"),
                  "Average effect (mean):",
                  min = -1, max = 1, value = 0.02, step = 0.01),
      helpText("The central effect size, around which sites vary."),
      sliderInput(paste0(prefix, "te_sd_normal"),
                  "Variation (SD):",
                  min = 0, max = 1, value = 0.01, step = 0.005),
      helpText("Controls how much sites differ from the average effect.")
    )
  )
}

ni_block_top <- function(prefix = "") {
  tagList(
    h4("Non-inferiority"),
    sliderInput(paste0(prefix, "ni_margin"),
                "Non-inferiority margin (0–0.3):",
                min = 0, max = 0.3, value = 0.10, step = 0.005),
    helpText("Sets how much worse the treatment can be compared to control 
             and still be considered acceptable. Smaller margins are stricter.")
  )
}

# Tab UI factory
tab_ui <- function(prefix, title, include_ni = FALSE) {
  tabPanel(
    title = title,
    value = prefix,
    fluidPage(
      titlePanel(if (include_ni) "Non-inferiority — Power / FP Visualizer"
                 else "Superiority — Power / FP Visualizer"),
      sidebarLayout(
        sidebarPanel(
          width = 5,
          h4("Core setup"),
          sliderInput(paste0(prefix, "alpha"),
                      "Alpha (one-sided alpha cutoff):",
                      min = 0, max = 0.5, value = 0.05, step = 0.025),
          helpText("This is the probability of incorrectly finding an effect when there isn’t one. 
                   For example, alpha = 0.05 means about a 5% risk of a false positive."),
          numericInput(paste0(prefix, "okay_sd"),
                       "Acceptable % error SD in power/false positive rate calculation:",
                       value = 1, min = 0.1, step = 1),
          helpText("This sets how precise you want the simulation results to be. 
          A smaller value means more accurate estimates of power or false positive rate, 
          but it will take longer to run. Larger values finish faster but are less precise."),
          tags$hr(),
          
          h4("Sample size range"),
          fluidRow(
            column(4, numericInput(paste0(prefix, "n_min"), "Minimum Sample Size:", value = 840, min = 2, step = 1)),
            column(4, numericInput(paste0(prefix, "n_max"), "Maximum Sample Size:", value = 1000, min = 2, step = 1)),
            column(4, numericInput(paste0(prefix, "n_by"),  "By:",    value = 10,  min = 1, step = 1))
          ),
          helpText("These settings control the range of sample sizes to test in the simulations. 
          For example, Min=840, Max=1000, By=10 will run simulations at 840, 850, 860, …, 1000."),
          tags$hr(),
          
          h4("Site counts"),
          textInput(paste0(prefix, "site_counts"),
                    "Site counts to simulate (comma separated ≥ 2):",
                    placeholder = "e.g., 2, 4, 6, 8, 10"),
          
          tags$hr(),
          site_mortality_block(prefix),
          
          sliderInput(paste0(prefix, "individual_sd"),
                      "Individual SD (0–1):",
                      min = 0, max = 1, value = 0.05, step = 0.005),
          helpText("How much patients within a site differ from the site's average mortality. 
          Higher values = more variation, lower values = patients are more similar."),
          tags$hr(),
          recruitment_block(prefix),
          tags$hr(),
          treatment_effect_block(prefix),
          if (include_ni) {
            tagList(tags$hr(), ni_block_top(prefix))
          } else tagList(),
          tags$hr(),
          uiOutput(paste0(prefix, "calcUI"))
        ),
        mainPanel(
          width = 7,
          h4("Status"),
          verbatimTextOutput(paste0(prefix, "status")),
          conditionalPanel(
            condition = sprintf("input['%scalc'] > 0", prefix),
            h4("Results"),
            plotOutput(paste0(prefix, "plot_power"), height = "340px"),
            plotOutput(paste0(prefix, "plot_fpr"),   height = "340px"),
            tags$hr()
          )
        )
      )
    )
  )
}

# ---------- UI ----------
ui <- navbarPage(
  title = "Visualizing the Impact of Site Heterogeneity on Power and Type I Error",
  id = "top_tabs",
  tab_ui(prefix = "sup_", title = "Superiority", include_ni = FALSE),
  tab_ui(prefix = "ni_",  title = "Non-inferiority", include_ni = TRUE)
)

# ---------- Server ----------
server <- function(input, output, session) {
  
  # ----- dynamic UI builders -----
  build_mortality_lists_ui <- function(prefix) {
    site_counts <- unique(na.omit(parse_int_list(nz(input[[paste0(prefix, "site_counts")]], na = ""))))
    site_counts <- site_counts[site_counts >= 2]
    if (length(site_counts) == 0) return(NULL)
    div(
      lapply(site_counts, function(s) {
        wellPanel(
          h5(sprintf("Site baseline mortalities — %d-site case (comma separated, values in [0,1])", s)),
          textInput(
            inputId   = paste0(prefix, "mortality_list_", s),
            label     = NULL,
            placeholder = paste(rep("0.10", s), collapse = ", ")
          )
        )
      })
    )
  }
  
  build_recruit_props_ui <- function(prefix) {
    site_counts <- unique(na.omit(parse_int_list(nz(input[[paste0(prefix, "site_counts")]], na = ""))))
    site_counts <- site_counts[site_counts >= 2]
    if (length(site_counts) == 0) return(NULL)
    div(
      lapply(site_counts, function(s) {
        wellPanel(
          h5(sprintf("Recruitment proportions — %d-site case (comma separated, sum to 1)", s)),
          textInput(
            inputId   = paste0(prefix, "recruit_props_", s),
            label     = NULL,
            placeholder = paste0(
              paste(rep(format(1/s, digits = 3), s), collapse = ", ")
            )
          )
        )
      })
    )
  }
  
  # render dynamic containers for both tabs
  for (prefix in c("sup_", "ni_")) {
    local({
      pfx <- prefix
      output[[paste0(pfx, "mortality_lists_container")]] <- renderUI(build_mortality_lists_ui(pfx))
      output[[paste0(pfx, "recruit_props_container")]]   <- renderUI(build_recruit_props_ui(pfx))
    })
  }
  
  # ----- shared readers/validators (prefix-aware) -----
  read_common <- function(prefix = "") {
    # Gather base inputs
    p <- list(
      alpha          = nz(input[[paste0(prefix, "alpha")]]),
      okay_sd        = nz(input[[paste0(prefix, "okay_sd")]]),
      
      n_min          = nz(input[[paste0(prefix, "n_min")]]),
      n_max          = nz(input[[paste0(prefix, "n_max")]]),
      n_by           = nz(input[[paste0(prefix, "n_by")]]),
      
      site_counts    = unique(na.omit(parse_int_list(nz(input[[paste0(prefix, "site_counts")]], na = "")))),
      
      mortality_mode = nz(input[[paste0(prefix, "mortality_mode")]], na = NA_character_),
      mortality_mean = nz(input[[paste0(prefix, "mortality_mean")]]),
      mortality_sd   = nz(input[[paste0(prefix, "mortality_sd")]]),
      mortality_range = nz(input[[paste0(prefix, "mortality_range")]]),
      
      individual_sd  = nz(input[[paste0(prefix, "individual_sd")]]),
      
      recruit_mode   = nz(input[[paste0(prefix, "recruit_mode")]], na = NA_character_),
      recruit_sd     = nz(input[[paste0(prefix, "recruit_sd")]]),
      
      te_mode            = nz(input[[paste0(prefix, "te_mode")]], na = NA_character_),
      te_const           = nz(input[[paste0(prefix, "te_const")]]),
      
      # new distinct inputs
      te_center_uniform  = nz(input[[paste0(prefix, "te_center_uniform")]]),
      te_width_uniform   = nz(input[[paste0(prefix, "te_width_uniform")]]),
      te_mean_normal     = nz(input[[paste0(prefix, "te_mean_normal")]]),
      te_sd_normal       = nz(input[[paste0(prefix, "te_sd_normal")]])
    )
    
    # Mortality lists by site count (only if per-site mode)
    if (identical(p$mortality_mode, "per_site") && length(p$site_counts) > 0) {
      ml <- setNames(vector("list", length(p$site_counts)), as.character(p$site_counts))
      for (s in p$site_counts) {
        ml[[as.character(s)]] <- parse_num_list(nz(input[[paste0(prefix, "mortality_list_", s)]], na = ""))
      }
      p$mortality_lists <- ml
    } else {
      p$mortality_lists <- list()
    }
    
    # Recruitment custom props by site count (if selected)
    if (identical(p$recruit_mode, "custom_props") && length(p$site_counts) > 0) {
      rp <- setNames(vector("list", length(p$site_counts)), as.character(p$site_counts))
      for (s in p$site_counts) {
        rp[[as.character(s)]] <- parse_num_list(nz(input[[paste0(prefix, "recruit_props_", s)]], na = ""))
      }
      p$recruit_props_lists <- rp
    } else {
      p$recruit_props_lists <- list()
    }
    
    p
  }
  
  validate_common <- function(p, label = "") {
    msgs <- c()
    
    # alpha, okay_sd
    if (is.na(p$alpha) || p$alpha <= 0 || p$alpha >= 1)
      msgs <- c(msgs, sprintf("%sOne-sided alpha must be (0,1).", label))
    if (is.na(p$okay_sd) || p$okay_sd < 0.1 || p$okay_sd > 10)
      msgs <- c(msgs, sprintf("%sAcceptable percent error SD must be in [0.1, 10].", label))
    
    # sample size range
    if (is.na(p$n_min) || p$n_min < 2 || p$n_min %% 1 != 0)
      msgs <- c(msgs, sprintf("%sMin N must be integer ≥ 2.", label))
    if (is.na(p$n_max) || p$n_max < 2 || p$n_max %% 1 != 0)
      msgs <- c(msgs, sprintf("%sMax N must be integer ≥ 2.", label))
    if (is.na(p$n_by)  || p$n_by  < 1 || p$n_by  %% 1 != 0)
      msgs <- c(msgs, sprintf("%sBy step must be integer ≥ 1.", label))
    if (length(msgs) == 0 && p$n_max <= p$n_min)
      msgs <- c(msgs, sprintf("%sMax N must be greater than Min N.", label))
    
    # site counts
    sc <- p$site_counts
    if (length(sc) == 0) {
      msgs <- c(msgs, sprintf("%sEnter at least one site count (≥ 2).", label))
    } else {
      if (any(is.na(sc)) || any(sc < 2) || any(sc %% 1 != 0))
        msgs <- c(msgs, sprintf("%sAll site counts must be integers ≥ 2.", label))
    }
    
    # mortality
    if (identical(p$mortality_mode, "per_site")) {
      for (s in sc) {
        v <- p$mortality_lists[[as.character(s)]] %||% numeric(0)
        if (length(v) != s)
          msgs <- c(msgs, sprintf("%sMortality list length for %d-site case (%d) must equal %d.", label, s, length(v), s))
        if (length(v) > 0 && (any(is.na(v)) || min(v) < 0 || max(v) > 1))
          msgs <- c(msgs, sprintf("%sMortality values for %d-site case must be numeric in [0,1].", label, s))
      }
    } else if (identical(p$mortality_mode, "mean_sd")) {
      if (is.na(p$mortality_mean) || p$mortality_mean < 0 || p$mortality_mean > 1)
        msgs <- c(msgs, sprintf("%sOverall mortality must be in [0,1].", label))
      if (is.na(p$mortality_sd) || p$mortality_sd < 0 || p$mortality_sd > 1)
        msgs <- c(msgs, sprintf("%sBetween-site SD must be in [0,1].", label))
    } else if (identical(p$mortality_mode, "uniform_range")) {
      r <- p$mortality_range
      if (length(r) != 2 || any(is.na(r)))
        msgs <- c(msgs, sprintf("%sUniform range: please set both Min and Max.", label))
      else {
        if (r[1] < 0 || r[1] > 1) msgs <- c(msgs, sprintf("%sUniform range: Min must be in [0,1].", label))
        if (r[2] < 0 || r[2] > 1) msgs <- c(msgs, sprintf("%sUniform range: Max must be in [0,1].", label))
        if (r[2] < r[1])         msgs <- c(msgs, sprintf("%sUniform range: Max must be ≥ Min.", label))
      }
    } else {
      msgs <- c(msgs, sprintf("%sChoose a mortality input mode.", label))
    }
    
    # individual SD
    if (is.na(p$individual_sd) || p$individual_sd < 0 || p$individual_sd > 1)
      msgs <- c(msgs, sprintf("%sIndividual SD must be in [0,1].", label))
    
    # recruitment
    if (identical(p$recruit_mode, "balanced")) {
      # ok
    } else if (identical(p$recruit_mode, "imbalance_sd")) {
      if (is.na(p$recruit_sd) || p$recruit_sd < 0 || p$recruit_sd > 1)
        msgs <- c(msgs, sprintf("%sRecruitment SD must be in [0,1].", label))
    } else if (identical(p$recruit_mode, "custom_props")) {
      for (s in sc) {
        r <- p$recruit_props_lists[[as.character(s)]] %||% numeric(0)
        if (length(r) != s)
          msgs <- c(msgs, sprintf("%sRecruitment vector length for %d-site case (%d) must equal %d.", label, s, length(r), s))
        if (length(r) > 0 && (any(is.na(r)) || min(r) < 0 || max(r) > 1))
          msgs <- c(msgs, sprintf("%sRecruitment proportions for %d-site case must be in [0,1].", label, s))
        if (length(r) > 0 && abs(sum(r) - 1) > 1e-6)
          msgs <- c(msgs, sprintf("%sRecruitment proportions for %d-site case must sum to 1 (currently %.6f).", label, s, sum(r)))
      }
    } else {
      msgs <- c(msgs, sprintf("%sChoose a recruitment mode.", label))
    }
    
    # treatment effect — validate mode-specific controls
    if (identical(p$te_mode, "constant")) {
      if (is.na(p$te_const) || p$te_const < -1 || p$te_const > 1)
        msgs <- c(msgs, sprintf("%sConstant effect must be in [-1,1].", label))
    } else if (identical(p$te_mode, "uniform")) {
      if (is.na(p$te_center_uniform) || p$te_center_uniform < -1 || p$te_center_uniform > 1)
        msgs <- c(msgs, sprintf("%sUniform center must be in [-1,1].", label))
      if (is.na(p$te_width_uniform) || p$te_width_uniform < 0 || p$te_width_uniform > 1)
        msgs <- c(msgs, sprintf("%sUniform width must be in [0,1].", label))
    } else if (identical(p$te_mode, "normal")) {
      if (is.na(p$te_mean_normal) || p$te_mean_normal < -1 || p$te_mean_normal > 1)
        msgs <- c(msgs, sprintf("%sNormal mean must be in [-1,1].", label))
      if (is.na(p$te_sd_normal) || p$te_sd_normal < 0 || p$te_sd_normal > 1)
        msgs <- c(msgs, sprintf("%sNormal SD must be in [0,1].", label))
    } else {
      msgs <- c(msgs, sprintf("%sChoose a treatment-effect mode.", label))
    }
    
    msgs
  }
  
  # ----- tab registration (wires a full tab by prefix) -----
  register_tab <- function(prefix, is_ni_tab = FALSE) {
    # params
    params <- reactive({
      base <- read_common(prefix)
      base$is_ni     <- is_ni_tab
      base$ni_margin <- if (is_ni_tab) nz(input[[paste0(prefix, "ni_margin")]]) else NA_real_
      
      # Derive mode-specific effective values for the simulator
      base$te_mean_eff   <- if (identical(base$te_mode, "uniform")) base$te_center_uniform
      else if (identical(base$te_mode, "normal"))  base$te_mean_normal
      else                                         NA_real_
      base$te_spread_eff <- if (identical(base$te_mode, "uniform")) base$te_width_uniform else NA_real_
      base$te_sd_eff     <- if (identical(base$te_mode, "normal"))  base$te_sd_normal     else NA_real_
      
      base
    })
    
    # validation
    validate_tab <- reactive({
      p <- params()
      msgs <- validate_common(p, label = "")
      if (isTRUE(p$is_ni)) {
        if (is.na(p$ni_margin) || p$ni_margin < 0 || p$ni_margin > 0.3)
          msgs <- c(msgs, "NI margin must be in [0, 0.3].")
      }
      msgs
    })
    
    # status
    output[[paste0(prefix, "status")]] <- renderText({
      msgs <- validate_tab()
      if (length(msgs) == 0) {
        "All inputs look valid. Click 'Run simulation' to proceed."
      } else {
        paste("Please address the following:", paste(paste0("- ", msgs), collapse = "\n"), sep = "\n")
      }
    })
    
    # calc button UI
    output[[paste0(prefix, "calcUI")]] <- renderUI({
      if (length(validate_tab()) == 0) {
        actionButton(paste0(prefix, "calc"), "Run simulation", class = "btn btn-primary")
      } else {
        tags$div(class = "text-muted", "Fix inputs to enable the Run simulation button.")
      }
    })
    
    # simulation
    sim_results <- eventReactive(input[[paste0(prefix, "calc")]], {
      p <- params()
      
      # derive grid
      Ns <- seq(from = as.integer(p$n_min), to = as.integer(p$n_max), by = as.integer(p$n_by))
      site_counts <- sort(unique(p$site_counts))
      # iteration scaling from acceptable SD
      baseline_iters <- 1000
      baseline_sd    <- 1.28
      iter_factor <- (baseline_sd / p$okay_sd)^2
      n_trials    <- max(100, round(baseline_iters * iter_factor))
      
      source("simulate_trial_rshiny.R")  # expects simulate_power_fp()
      
      total <- length(Ns) * length(site_counts)
      progressed <- 0
      step <- if (total > 0) 1/total else 1
      
      rows <- list()
      idx <- 0L
      
      withProgress(message = "Running simulations…", value = 0, {
        for (s in site_counts) {
          # decide mortality mode per s
          mort_list_s <- NULL
          mmode <- as.character(p$mortality_mode)
          m_mean <- NA_real_
          m_sd   <- NA_real_
          
          if (identical(mmode, "per_site")) {
            mort_list_s <- p$mortality_lists[[as.character(s)]]
          } else if (identical(mmode, "mean_sd")) {
            m_mean <- p$mortality_mean
            m_sd   <- p$mortality_sd
          } else if (identical(mmode, "uniform_range")) {
            # convert to per-site list with evenly spaced values from range slider
            mmode <- "per_site"
            r <- p$mortality_range
            mort_list_s <- seq(r[1], r[2], length.out = s)
          }
          
          # recruitment props for custom mode
          recruit_props_s <- if (identical(p$recruit_mode, "custom_props")) {
            p$recruit_props_lists[[as.character(s)]]
          } else numeric(0)
          
          for (n in Ns) {
            #print(p$te_mean_eff)
            sim <- simulate_power_fp(
              sample_size    = n,
              n_sites        = s,
              
              mortality_mode = mmode,
              mortality_list = if (identical(mmode, "per_site")) mort_list_s else NULL,
              mortality_mean = if (identical(mmode, "mean_sd")) m_mean else NA_real_,
              mortality_sd   = if (identical(mmode, "mean_sd")) m_sd   else NA_real_,
              
              individual_sd  = p$individual_sd,
              
              recruit_mode   = as.character(p$recruit_mode),
              recruit_sd     = if (identical(p$recruit_mode, "imbalance_sd")) p$recruit_sd else NA_real_,
              recruit_props  = recruit_props_s,
              
              te_mode        = as.character(p$te_mode),
              te_const       = p$te_const,
              te_mean        = p$te_mean_eff,
              te_spread      = p$te_spread_eff,
              te_sd          = p$te_sd_eff,
              
              is_ni          = isTRUE(p$is_ni),
              ni_margin      = if (isTRUE(p$is_ni)) p$ni_margin else NA_real_,
              
              n_trials       = n_trials,
              alpha          = p$alpha
            )
            
            idx <- idx + 1L
            rows[[idx]] <- data.frame(
              sample_size = as.integer(n),
              n_sites     = as.integer(s),
              power       = as.numeric(sim$power),  # expected in %
              fpr         = as.numeric(sim$fpr),    # expected in %
              stringsAsFactors = FALSE
            )
            
            progressed <- min(1, progressed + step)
            incProgress(step, detail = sprintf("Sites = %d, Sample Size = %d  (%d/%d)", s, n, idx, total))
          }
        }
      })
      
      df <- if (length(rows)) do.call(rbind, rows) else
        data.frame(sample_size = integer(0), n_sites = integer(0), power = numeric(0), fpr = numeric(0))
      
      list(df = df, n_trials = n_trials)
    })
    
    # ----- plots -----
    output[[paste0(prefix, "plot_power")]] <- renderPlot({
      req(input[[paste0(prefix, "calc")]] > 0)
      res <- sim_results()
      df  <- res$df
      validate(need(nrow(df) > 0, "No results to display."))
      
      ggplot(df, aes(x = sample_size, y = power, group = factor(n_sites), color = factor(n_sites))) +
        geom_point(alpha = 0.7) +
        ggplot2::geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = TRUE,
          span = 0.75,
          method.args = list(degree = 1, family = "gaussian"),
          linewidth = 1,
          alpha = 0.35,
          na.rm = TRUE
        ) +
        labs(x = "Sample size (total, both arms)",
             y = "Power (%)",
             color = "Sites",
             title = sprintf("Power vs Sample Size (Iterations per point: %d)", res$n_trials)) +
        theme_minimal(base_size = 12)
    })
    
    output[[paste0(prefix, "plot_fpr")]] <- renderPlot({
      req(input[[paste0(prefix, "calc")]] > 0)
      res <- sim_results()
      df  <- res$df
      validate(need(nrow(df) > 0, "No results to display."))
      
      ggplot(df, aes(x = sample_size, y = fpr, group = factor(n_sites), color = factor(n_sites))) +
        geom_point(alpha = 0.7) +
        ggplot2::geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = TRUE,
          span = 0.75,
          method.args = list(degree = 1, family = "gaussian"),
          linewidth = 1,
          alpha = 0.35,
          na.rm = TRUE
        ) +
        labs(x = "Sample size (total, both arms)",
             y = "False positive rate (%)",
             color = "Sites",
             title = sprintf("False Positive Rate vs Sample Size (n_trials per point: %d)", res$n_trials)) +
        theme_minimal(base_size = 12)
    })
  }
  
  # Register both tabs
  register_tab(prefix = "sup_", is_ni_tab = FALSE)
  register_tab(prefix = "ni_",  is_ni_tab = TRUE)
}

shinyApp(ui, server)
