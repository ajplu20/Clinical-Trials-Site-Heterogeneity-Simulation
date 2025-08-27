library(ggplot2)
library(viridis)       # For colorblind-friendly palettes
library(dplyr)
library(tidyr)
library(patchwork)
library(glue)
library(scales)


my_colors <- c("#451077FF","#365C8DFF","#4AC16DFF", "#FD9567FF")

plot_test_power_by_site <- function(results_df, individual_sd, treatment_sd, recruitment_sd, title_num) {
  # Define mapping from test variable names to human-readable labels
  test_labels <- c(
    naive_t_test = "Naive t-test",
    fix_effect_model_test = "Fixed effect model",
    mixed_effect_model_test = "Mixed effect model",
    mantel_haenszel_test = "Mantel-Haenszel test"
  )
  
  # Tests to extract
  test_names <- names(test_labels)
  
  # Prepare summary dataframe: one row per (site, test)
  summary_df <- lapply(seq_len(nrow(results_df)), function(i) {
    site_num <- results_df$num_sites[i]
    sim_data <- results_df$simulation_data[[i]]
    
    means <- sapply(test_names, function(test) {
      mean(sim_data[[test]], na.rm = TRUE)
    })
    
    data.frame(num_sites = site_num, t(means))
  }) %>% bind_rows()
  
  # Convert to long format for ggplot and relabel test types
  plot_data <- summary_df %>%
    pivot_longer(
      cols = all_of(test_names),
      names_to = "test_type",
      values_to = "mean_result"
    ) %>%
    mutate(test_type = recode(test_type, !!!test_labels))
  
  # Plot title with parameter values
  plot_title = case_when(
    title_num == 1 ~ "Variable Treatment Effect per Individual Patient",
    title_num == 2 ~ "Variable Treatment Effect per Study Site",
    title_num == 3 ~ "Variable Number of Patients Enrolled at Each Study Site",
    title_num == 4 ~ "Variable Treatment Effect and Number of Patients Enrolled at Each Study Site",
    TRUE           ~ NA_character_    # fallback if title_num is missing or out of range
  )

  # Final ggplot
  return(
    ggplot(plot_data, aes(x = num_sites, y = mean_result, color = test_type)) +
      geom_line(alpha = 0.5, size = 1) +
      geom_point(size = 2) +
      labs(
        title = plot_title,
        x = "Number of Sites",
        y = "Power",
        color = "Analysis Method"
      ) + 
      ylim(0.5, 1) +
      scale_color_manual(values = my_colors) +
      theme_minimal() +
      theme(text = element_text(size = 12)) +
      geom_hline(yintercept = 0.9, linetype = "dotted", color = "black")
  )
}



plot_power_thresholds_by_sd <- function(df, col) {
  # Reshape to long format
  plot_data <- df %>%
    pivot_longer(
      cols = -{{col}},
      names_to = "test_type",
      values_to = "num_sites"
    ) %>%
    mutate(
      test_type = recode(test_type,
                         naive_t_test = "Naive t-test",
                         fix_effect_model_test = "Fixed effect model",
                         mixed_effect_model_test = "Mixed effect model",
                         mantel_haenszel_test = "Mantel-Haenszel test"
      )
    )
  
  # Plot
  return(
    ggplot(plot_data, aes(x = .data[[col]], y = num_sites, color = test_type)) +
      geom_line(alpha = 0.5, size = 1) +
      geom_point(size = 2) +
      labs(
        title = "",
        x = "Standard Deviation",
        y = "Number of Sites Required to Achieve 90% Power",
        color = "Analysis Method"
      ) +
      scale_color_manual(values = my_colors) +
      scale_y_continuous(breaks = seq(0, 100, by = 5), limits = c(1, 100)) +
      theme_minimal() +
      theme(text = element_text(size = 12))
  )
}



plot_power_thresholds_by_mid_50_percent <- function(df, col) {
  # Reshape to long format
  plot_data <- df %>%
    pivot_longer(
      cols = -{{col}},
      names_to = "test_type",
      values_to = "num_sites"
    ) %>%
    mutate(
      test_type = recode(test_type,
                         naive_t_test = "Naive t-test",
                         fix_effect_model_test = "Fixed effect model",
                         mixed_effect_model_test = "Mixed effect model",
                         mantel_haenszel_test = "Mantel-Haenszel test"
      ),
      # Calculate 25th and 75th percentiles of mortality based on individual_sd
      lower = qnorm(0.25, mean = 0, sd = .data[[col]]),
      upper = qnorm(0.75, mean = 0, sd = .data[[col]]),
      mortality_range = paste0("-", percent(lower, accuracy = 1), " to ", percent(upper, accuracy = 1))
    )
  
  # To keep x-axis ordered nicely
  plot_data$mortality_range <- factor(plot_data$mortality_range, levels = unique(plot_data$mortality_range))
  
  # Plot
  return(
    ggplot(plot_data, aes(x = mortality_range, y = num_sites, color = test_type, group = test_type)) +
      geom_line(alpha = 0.5, size = 1) +
      geom_point(size = 2) +
      labs(
        title = "",
          x = "Middle 50% Range",
        y = "Number of Sites Required to Achieve 90% Power",
        color = "Analysis Method"
      ) +
      scale_color_manual(values = my_colors) +
      scale_y_continuous(breaks = seq(0, 100, by = 5), limits = c(1, 100)) +
      theme_minimal() +
      theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))
  )
}


#0.1,0,0, only varying the individual sd
result100_df <- readRDS("individualsd=0.1_treatmentsd = 0.0_recruitmentsd=0.0.rds")
p100Left <- plot_test_power_by_site(result100_df, 0.1,0,0,1)
p100Left

power_df_individual <- readRDS("individualsd_treatmentsd=0_recruitmentsd=0.rds")
p100Right_SD <- plot_power_thresholds_by_sd(power_df_individual, "individual_sd")
p100Right_SD

p100Right_50_percent <- plot_power_thresholds_by_mid_50_percent(power_df_individual, "individual_sd")
p100Right_50_percent



#0.05, 0.1,0, varying the treatment sd,
result100_df <- readRDS("individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.0.rds")
p510Left <- plot_test_power_by_site(result100_df, 0.1,0,0,2)
p510Left

power_df_individual <- readRDS("treatmentsd_individualsd=0.05_recruitmentsd=0.rds")
power_df_individual <- power_df_individual[-nrow(power_df_individual), ]

p510Right_SD <- plot_power_thresholds_by_sd(power_df_individual, "treatment_sd")
p510Right_SD

p510Right_50_percent <- plot_power_thresholds_by_mid_50_percent(power_df_individual, "treatment_sd")
p510Right_50_percent


combined <- (p100Left + p100Right_50_percent + p510Left + p510Right_50_percent) +
  plot_layout(nrow = 1, guides = "collect") & 
  theme(legend.position = "bottom")

print(combined)

#0.05 ,0,0.1, varying the recruitment sd, 

result100_df <- readRDS("individualsd=0.05_treatmentsd = 0.0_recruitmentsd=0.1.rds")
p501Left <- plot_test_power_by_site(result100_df, 0.1,0,0,3)
p501Left

power_df_individual <- readRDS("recruitmentsd_individualsd=0.05_treatmentsd=0.rds")

p501Right_SD <- plot_power_thresholds_by_sd(power_df_individual, "recruitment_sd")
p501Right_SD

p501Right_50_percent <- plot_power_thresholds_by_mid_50_percent(power_df_individual, "recruitment_sd")
p501Right_50_percent




#0.05,0.1,0.1, varying the recruitment sd, with 0.1 for the other 2.

result100_df <- readRDS("individualsd=0.05_treatmentsd = 0.1_recruitmentsd=0.1.rds")
p511Left <- plot_test_power_by_site(result100_df, 0.1,0,0,2)
p511Left

power_df_individual <- readRDS("treatmentsd_individualsd=0.1_recruitmentsd=0.1.rds")
p511Right_SD <- plot_power_thresholds_by_sd(power_df_individual, "treatment_sd")
p511Right_SD

p511Right_50_percent <- plot_power_thresholds_by_mid_50_percent(power_df_individual, "treatment_sd")
p511Right_50_percent





combined <- (
  p100Left + p100Right_50_percent + 
    p510Left + p510Right_50_percent +  
    p501Left + p501Right_50_percent +  
    p511Left + p511Right_50_percent
) +
  plot_layout(ncol = 2, nrow = 4, guides = "collect") & 
  theme(legend.position = "bottom")

print(combined)

ggsave("combined_plot.png", combined, width = 12, height = 24, dpi = 300)
