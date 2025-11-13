# Save as app.R and run with: shiny::runApp()

library(shiny)
library(ggplot2)

# ========= easy-to-tune hyperparameters =========
N_TRIALS_DEFAULT   <- 1000   # iterations per (N, sites) point
GRID_STEPS_DEFAULT <- 10     # <- NEW: how many steps between n_min and n_max

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
site_mortality_block <- function(prefix = "", show_header = TRUE) {
  tagList(
    if (show_header) h4("Site baseline mortality (control group)"),
    selectInput(
      paste0(prefix, "mortality_mode"), "How site mean mortality differ:",
      choices = c("Manually enter the mean for each site" = "per_site",
                  "Evenly spaced across a specified range" = "uniform_range",
                  "Randomly drawn from a Normal(mean, SD) distribution" = "mean_sd"),
      selected = "per_site"
    ),
    
    # mean/sd mode
    conditionalPanel(
      condition = sprintf("input['%smortality_mode'] == 'mean_sd'", prefix),
      sliderInput(paste0(prefix, "mortality_mean"),
                  "mean:",
                  min = 0, max = 1, value = 0.10, step = 0.005),
      sliderInput(paste0(prefix, "mortality_sd"),
                  "SD:",
                  min = 0, max = 1, value = 0.03, step = 0.005)
    ),
    
    # per-site mode: dynamic rows for each site count
    conditionalPanel(
      condition = sprintf("input['%smortality_mode'] == 'per_site'", prefix),
      uiOutput(paste0(prefix, "mortality_lists_container")),
      helpText("Enter the mean control group mortality for each site directly.")
    ),
    
    # uniform range mode: range slider min/max
    conditionalPanel(
      condition = sprintf("input['%smortality_mode'] == 'uniform_range'", prefix),
      sliderInput(paste0(prefix, "mortality_range"),
                  "Mortality range (min to max):",
                  min = 0, max = 1, value = c(0.05, 0.15), step = 0.005),
      helpText("Each site's mean control group mortality will be spread evenly between the minimum and maximum mortality values.")
    )
  )
}

recruitment_block <- function(prefix = "") {
  tagList(
    h4("Sample size balance across sites"),
    selectInput(
      paste0(prefix, "recruit_mode"), "How to distribute patients across sites:",
      choices = c("Randomized imbalance around equal recruitment:" = "imbalance_sd",
                  "Manually enter proportion of patients from each site" = "custom_props"),
      selected = "imbalance_sd"
    ),
    
    # imbalance_sd mode
    conditionalPanel(
      condition = sprintf("input['%srecruit_mode'] == 'imbalance_sd'", prefix),
      sliderInput(paste0(prefix, "recruit_sd"),
                  "Sampling imbalance (SD):", 
                  min = 0, max = 1, value = 0.05, step = 0.005),
      helpText("0 is equal number of patients recruited from each site. Larger means recruitment is more imbalanced across sites (some sites have more patients, some have less).")
    ),
    
    # custom proportions mode
    conditionalPanel(
      condition = sprintf("input['%srecruit_mode'] == 'custom_props'", prefix),
      uiOutput(paste0(prefix, "recruit_props_container")),
      helpText("Enter the exact recruitment proportion for each site (must sum to 1).")
    )
  )
}

treatment_effect_block <- function(prefix = "", show_header = TRUE) {
  tagList(
    h4("Treatment effect at each site"),
    helpText("site control group mortality - site treatment group mortality"),
    selectInput(
      paste0(prefix, "te_mode"), "Input mode:",
      choices = c("Constant (same at all sites)" = "constant",
                  "Uniform (drawn randomly within a range)" = "uniform",
                  "Normal (drawn randomly from N(mean, SD)" = "normal"),
      selected = "constant"
    ),
    
    # constant effect mode
    conditionalPanel(
      condition = sprintf("input['%ste_mode'] == 'constant'", prefix),
      sliderInput(paste0(prefix, "te_const"),
                  "Average treatment effect:",
                  min = -1, max = 1, value = 0.10, step = 0.01)
    ),
    
    # uniform effect mode — two-point range slider
    conditionalPanel(
      condition = sprintf("input['%ste_mode'] == 'uniform'", prefix),
      sliderInput(paste0(prefix, "te_range_uniform"),
                  "min and max of the range:",
                  min = -1, max = 1, value = c(-0.10, 0.30), step = 0.01)
    ),
    
    # normal effect mode
    conditionalPanel(
      condition = sprintf("input['%ste_mode'] == 'normal'", prefix),
      sliderInput(paste0(prefix, "te_mean_normal"),
                  "mean:",
                  min = -1, max = 1, value = 0.02, step = 0.01),
      sliderInput(paste0(prefix, "te_sd_normal"),
                  "SD:",
                  min = 0, max = 1, value = 0.01, step = 0.005)
    )
  )
}

ni_block_top <- function(prefix = "") {
  tagList(
    h4("Non-inferiority margin"),
    helpText("maximum acceptable % treatment is worse than control (maximum treatment mortality - control mortality value that is still considered not too bad)"),
    sliderInput(paste0(prefix, "ni_margin"),
                "select value in range (0–0.3):",
                min = 0, max = 0.3, value = 0.10, step = 0.005)
  )
}

# --- NEW: single heterogeneity section (mortality -> individual SD -> treatment)
heterogeneity_block <- function(prefix = "") {
  tagList(
    h4("Control group mortality"),
    site_mortality_block(prefix, show_header = FALSE),
    sliderInput(
      paste0(prefix, "individual_sd"),
      "How much patient mortality varies within each site",
      min = 0, max = 1, value = 0.05, step = 0.005
    ),
    treatment_effect_block(prefix, show_header = FALSE)
  )
}

# Tab UI factory (ORDERED)
tab_ui <- function(prefix, title, include_ni = FALSE) {
  tabPanel(
    title = title,
    value = prefix,
    fluidPage(
      titlePanel(if (include_ni) "Non-inferiority — Power / Type 1 error visualizer"
                 else "Superiority — Power / Type 1 error visualizer"),
      sidebarLayout(
        sidebarPanel(
          width = 5,
          
          # 1) Alpha
          h4("Core setup"),
          selectInput(paste0(prefix, "alpha"),
                      "Alpha (one-sided alpha cutoff):",
                      choices = c("0.05" = 0.05, "0.025" = 0.025, "0.0125" = 0.0125),
                      selected = 0.05),
          helpText("Probability of incorrectly finding an effect when there isn’t one."),
          tags$hr(),
          
          # 2) Site counts
          h4("Number of site"),
          textInput(paste0(prefix, "site_counts"),
                    "Specify how many sites the trial can include:",
                    placeholder = "e.g., 2, 4, 6, 8, 10"),
          tags$hr(),
          
          # 3) Sample size range
          h4("Sample size range"),
          fluidRow(
            column(6, numericInput(paste0(prefix, "n_min"), "Minimum Sample Size:", value = 840, min = 2, step = 1)),
            column(6, numericInput(paste0(prefix, "n_max"), "Maximum Sample Size:", value = 1000, min = 2, step = 1))
          ),
          helpText("Calculate power and type 1 error from Min to Max sample size."),
          tags$hr(),
          
          # 4) Recruitment
          recruitment_block(prefix),
          tags$hr(),
          
          # 5) Heterogeneity (mortality -> individual SD -> treatment)
          heterogeneity_block(prefix),
          
          # NI margin at the bottom (only on NI tab)
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
  selected = "instructions",   # <- NEW: open on the Instructions tab by default
  tabPanel(
    title = "Instructions",
    value = "instructions",
    fluidPage(
      titlePanel("Instructions"),
      fluidRow(
        column(
          width = 12,
          h3("How to use this app"),
          tags$ol(
            tags$li(HTML("<b>Choose the trial type:</b> use the <i>Superiority</i> or <i>Non-inferiority</i> tab at the top.<br/>
                          • <i>Superiority trials</i>: you want to show the new treatment works better.<br/>
                          • <i>Non-inferiority trials</i>: you want to show the new treatment is not unacceptably worse than control (for example, if it is cheaper or easier to give).")),
            tags$li(HTML("<b>Set the significance level (alpha):</b> in <i>Core setup ▸ Alpha</i>, choose how strict you want the test to be. A smaller value means it is harder to declare success, but you are less likely to make a false claim that the treatment works.")),
            tags$li(HTML("<b>Enter the number of sites:</b> in <i>Number of site</i>, type one or more values, such as <code>2, 5, 10</code>. Each value represents a different design: for example, a 2-site trial versus a 10-site trial.")),
            tags$li(HTML("<b>Set the total sample size range:</b> choose the <i>Minimum</i> and <i>Maximum Sample Size</i> (total number of patients across both arms). The app will automatically try several sample sizes between these two values.")),
            tags$li(HTML("<b>Describe how patients are spread across sites:</b> under <i>Sample size balance across sites</i> you can:<br/>
                          • Let the app create a mild random imbalance around equal recruitment (more realistic).<br/>
                          • Or manually type the percentage of patients that will come from each site (these percentages must add up to 100%).")),
            tags$li(HTML("<b>Describe how sick patients are at each site (baseline mortality):</b> under <i>Control group mortality</i>, choose one option:<br/>
                          • <i>Manually enter the mean for each site</i>: you type a value between 0 and 1 for each site (for example, 0.10 means 10% mortality in the control arm).<br/>
                          • <i>Evenly spaced across a specified range</i>: you give a minimum and maximum mortality, and the app spreads sites evenly across this range.<br/>
                          • <i>Randomly drawn from a Normal(mean, SD)</i>: you provide an overall average mortality and a spread (SD); sites are drawn around that average.")),
            tags$li(HTML("<b>Set how much outcomes vary within each site:</b> use <i>How much patient mortality varies within each site</i> to describe how different patients at the same site can be. Bigger values mean more variation from patient to patient.")),
            tags$li(HTML("<b>Describe the treatment effect across sites:</b> under <i>Treatment effect at each site</i>, choose:<br/>
                          • <i>Constant</i>: the treatment has the same benefit at every site.<br/>
                          • <i>Uniform</i>: the benefit can vary from site to site within a range you specify.<br/>
                          • <i>Normal</i>: the benefit varies around an average value with some spread (SD).<br/>
                          Here, the treatment effect is defined as <i>control mortality − treatment mortality</i>, so a positive value means the treatment reduces mortality, and a negative value means treatment increases mortality.")),
            tags$li(HTML("<b>For non-inferiority trials only:</b> set the <i>Non-inferiority margin</i>.</br>
                          This margin is the largest extra risk you are willing to accept for the new treatment compared with control. If the new treatment is no worse than this margin, it can be considered non-inferior.")),
            tags$li(HTML("<b>Run the simulation:</b> check the <i>Status</i> box on the right. If there are any input problems (for example, percentages not adding to 100%), they will be listed there. Once everything is valid, the <i>Run simulation</i> button will appear. Click it to start."))
          ),
          tags$hr(),
          h3("What the graphs show"),
          tags$ul(
            tags$li(HTML("<b>Power vs Sample Size</b>: shows the chance that the trial will correctly detect a treatment effect, for each combination of sample size and number of sites. Each line corresponds to a different number of sites.")),
            tags$li(HTML("<b>Type 1 Error vs Sample Size</b>: shows how often the trial would <i>incorrectly</i> conclude that the treatment works (or is non-inferior) when in fact it does not.")),
            tags$li("By comparing the lines for different numbers of sites, you can see how adding more sites helps stabilise power and false positive rates.")
          ),
          tags$hr(),
          h3("How the simulation works (in simple terms)"),
          tags$p("For each design you specify (a choice of total sample size and number of sites), the app repeats the following steps many times:"),
          tags$ol(
            tags$li("It randomly assigns patients to treatment or control in a 1:1 ratio, without stratifying by site."),
            tags$li("It uses your inputs to decide how many patients are recruited at each site, how sick they are in the control arm, and how much benefit the treatment has at that site."),
            tags$li("It simulates who lives or dies at each site, under treatment and control."),
            tags$li("It analyses the trial and records whether the final conclusion was that the treatment works (or is non-inferior).")
          ),
          tags$p("By repeating this process many times (for example, 1,000 simulations per design), the app estimates:"),
          tags$ul(
            tags$li("the percentage of simulated trials that correctly conclude the treatment works (power), and"),
            tags$li("the percentage that incorrectly conclude the treatment works when it actually does not (false positive rate).")
          ),
          tags$p("These summaries help you understand how site-to-site differences and unbalanced recruitment might affect the reliability of your trial results."),
          tags$hr(),
          h3("Input rules and practical tips"),
          tags$ul(
            tags$li(HTML("<b>Site counts:</b> must be whole numbers 2 or larger, separated by commas (for example, <code>2, 4, 8</code>).")),
            tags$li(HTML("<b>Percentages and probabilities:</b> values such as mortality, proportions, and within-site variation must be between 0 and 1 (for example, 0.10 means 10%).")),
            tags$li(HTML("<b>Custom recruitment patterns:</b> when you enter recruitment proportions by site, make sure they add up to 1 (or 100%). The app will flag this if it is not satisfied.")),
            tags$li(HTML("<b>Sample size range:</b> choose a realistic minimum and maximum based on what is feasible in practice. The app will automatically pick several points in between to draw the lines on the graphs.")),
            tags$li(HTML("<b>Design insight:</b> in general, increasing the number of sites helps the test better, with higher power and lower type 1 error. Simply increasing the total number of patients does not always fix problems caused by large differences in treatment effect between sites."))
          ),
          tags$hr(),
          h3("Feedback"),
          tags$p(HTML("If you notice any issues or have suggestions, please contact <a href='mailto:mdcmy@nus.edu.sg'>Mo Yin</a>.")),
          h3("Latest update"),
          tags$p("12 November 2025 (SGT)")
        )
      )
      
      
    )
  ),
  tab_ui(prefix = "sup_", title = "Superiority trials", include_ni = FALSE),
  tab_ui(prefix = "ni_",  title = "Non-inferiority trials", include_ni = TRUE)
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
      alpha          = as.numeric(nz(input[[paste0(prefix, "alpha")]])),
      
      n_min          = nz(input[[paste0(prefix, "n_min")]]),
      n_max          = nz(input[[paste0(prefix, "n_max")]]),
      
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
      
      te_range_uniform   = nz(input[[paste0(prefix, "te_range_uniform")]]),
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
    
    # alpha
    if (is.na(p$alpha) || !(p$alpha %in% c(0.05, 0.025, 0.0125)))
      msgs <- c(msgs, sprintf("%sOne-sided alpha must be one of 0.05, 0.025, 0.0125.", label))
    
    # sample size range
    if (is.na(p$n_min) || p$n_min < 2 || p$n_min %% 1 != 0)
      msgs <- c(msgs, sprintf("%sMin N must be integer ≥ 2.", label))
    if (is.na(p$n_max) || p$n_max < 2 || p$n_max %% 1 != 0)
      msgs <- c(msgs, sprintf("%sMax N must be integer ≥ 2.", label))
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
    if (identical(p$recruit_mode, "imbalance_sd")) {
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
    
    # treatment effect
    if (identical(p$te_mode, "constant")) {
      if (is.na(p$te_const) || p$te_const < -1 || p$te_const > 1)
        msgs <- c(msgs, sprintf("%sConstant effect must be in [-1,1].", label))
    } else if (identical(p$te_mode, "uniform")) {
      r <- p$te_range_uniform
      if (length(r) != 2 || any(is.na(r)))
        msgs <- c(msgs, sprintf("%sUniform effect: please set both Min and Max.", label))
      else {
        if (r[1] < -1 || r[1] > 1) msgs <- c(msgs, sprintf("%sUniform effect: Min must be in [-1,1].", label))
        if (r[2] < -1 || r[2] > 1) msgs <- c(msgs, sprintf("%sUniform effect: Max must be in [-1,1].", label))
        if (r[2] < r[1])          msgs <- c(msgs, sprintf("%sUniform effect: Max must be ≥ Min.", label))
      }
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
  
  # ----- tab registration -----
  register_tab <- function(prefix, is_ni_tab = FALSE) {
    # params
    params <- reactive({
      base <- read_common(prefix)
      base$is_ni     <- is_ni_tab
      base$ni_margin <- if (is_ni_tab) nz(input[[paste0(prefix, "ni_margin")]]) else NA_real_
      
      # Derive mode-specific effective values for the simulator
      if (identical(base$te_mode, "uniform")) {
        rng <- base$te_range_uniform
        base$te_mean_eff   <- mean(rng)
        base$te_spread_eff <- diff(rng)
      } else if (identical(base$te_mode, "normal")) {
        base$te_mean_eff   <- base$te_mean_normal
        base$te_sd_eff     <- base$te_sd_normal
      } else {
        base$te_mean_eff   <- NA_real_
        base$te_spread_eff <- NA_real_
        base$te_sd_eff     <- NA_real_
      }
      
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
      
      # derive grid (auto step)
      step_auto <- max(1L, floor((as.integer(p$n_max) - as.integer(p$n_min)) / GRID_STEPS_DEFAULT))
      Ns <- seq(from = as.integer(p$n_min), to = as.integer(p$n_max), by = step_auto)
      if (tail(Ns, 1) != as.integer(p$n_max)) {
        Ns <- unique(c(Ns, as.integer(p$n_max)))
      }
      
      site_counts <- sort(unique(p$site_counts))
      n_trials    <- N_TRIALS_DEFAULT
      
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
            mmode <- "per_site"
            r <- p$mortality_range
            mort_list_s <- seq(r[1], r[2], length.out = s)
          }
          
          # recruitment props for custom mode
          recruit_props_s <- if (identical(p$recruit_mode, "custom_props")) {
            p$recruit_props_lists[[as.character(s)]]
          } else numeric(0)
          
          for (n in Ns) {
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
          span = 1.2,
          method.args = list(degree = 2, family = "gaussian"),
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
          span = 2.5,
          method.args = list(degree = 2, family = "gaussian"),
          linewidth = 1,
          alpha = 0.35,
          na.rm = TRUE
        ) +
        labs(x = "Sample size (total, both arms)",
             y = "Type 1 error (%)",
             color = "Sites",
             title = sprintf("Type 1 Error vs Sample Size (Iterations per point: %d)", res$n_trials)) +
        theme_minimal(base_size = 12)
    })
  }
  
  # Register both tabs
  register_tab(prefix = "sup_", is_ni_tab = FALSE)
  register_tab(prefix = "ni_",  is_ni_tab = TRUE)
}

shinyApp(ui, server)
