setwd("./sample_environments/environments/environment_final")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(cowplot)
library(gridExtra)
library(grid)
library(tidyr)
# Set theme for consistency
theme_set(theme_classic())

# ----------------------------
# 1. Data Preparation
# ----------------------------

# Load the master CSV file
file_path <- "./sample_environments/environments/environment_final/master_genome_coverage_all_final_hmp_removed.csv"
data <- read.csv(file_path)

# Count total number of samples per environment for all calculations
total_samples_per_env <- data %>%
  group_by(Environment) %>%
  summarise(Total_Samples = n_distinct(SampleID), .groups = 'drop')


# Filter data based on your threshold for presence
filtered_data <- data %>%
  filter(Breadth >= 0.6, Coverage >= 1)

# Calculate collective prevalence for each environment
collective_prevalence <- filtered_data %>%
  group_by(Environment) %>%
  summarise(Total_Present = n_distinct(SampleID), .groups = 'drop') %>%
  left_join(total_samples_per_env, by = "Environment") %>%
  mutate(Prevalence = (Total_Present / Total_Samples) * 100)

# Calculate overall total present and total samples
overall_total_present <- filtered_data %>%
  summarise(Total_Present = n_distinct(SampleID)) %>%
  pull(Total_Present)

overall_total_samples <- sum(total_samples_per_env$Total_Samples)

# Calculate overall prevalence
overall_prevalence <- 100 * (overall_total_present / overall_total_samples)

# Append overall data to collective prevalence
overall_data <- tibble(
  Environment = "overall",
  Total_Present = overall_total_present,
  Total_Samples = overall_total_samples,
  Prevalence = overall_prevalence
)

collective_prevalence <- bind_rows(collective_prevalence, overall_data)

# Calculate prevalence for each MAG within each environment
prevalence_data <- filtered_data %>%
  group_by(Environment, MAG) %>%
  summarise(Count = n_distinct(SampleID), .groups = 'drop') %>%
  left_join(total_samples_per_env, by = "Environment") %>%
  mutate(Prevalence = (Count / Total_Samples) * 100)

# ----------------------------
# 2. Visualization Components
# ----------------------------

# Function to create individual pie charts
create_pie_chart <- function(env, prevalence, total_samples, color) {
  df <- data.frame(
    Category = c("Prevalence", "Remaining"),
    Value = c(prevalence, 100 - prevalence)
  )
  
  ggplot(df, aes(x = "", y = Value, fill = Category)) +
    geom_col(width = 1, color = "black") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Prevalence" = color, "Remaining" = "darkgrey")) +
    labs(
      title = env,
      subtitle = paste0("Prevalence: ", formatC(prevalence, format = "f", digits = 2), "%\nn = ", total_samples)
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
}

# Generate pie charts including the overall data
pie_colors <- brewer.pal(n = length(unique(collective_prevalence$Environment)), name = "Set2")
pie_charts <- mapply(create_pie_chart, env = collective_prevalence$Environment, prevalence = collective_prevalence$Prevalence, total_samples = collective_prevalence$Total_Samples, color = pie_colors, SIMPLIFY = FALSE)
pie_chart_grid <- plot_grid(plotlist = pie_charts, ncol = 4)
print(pie_chart_grid)
output_path <- "/mnt/nrdstor/yinlab/jakresi/sample_environments/environments/figs/"
ggsave(file.path(output_path, "prevalence_analysis_completed_piechard_only_hmp_removed.pdf"), plot = pie_chart_grid, width = 11, height = 5, dpi = 300)

# Bar Plot
bar_plot <- ggplot(prevalence_data, aes(x = MAG, y = Prevalence, fill = Environment)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "MAG", y = "Prevalence (%)", title = "Prevalence of Each MAG in Each Environment") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(bar_plot)
output_path <- "./sample_environments/environments/figs/"
ggsave(file.path(output_path, "prevalence_analysis_completed_barplot_only.pdf"), plot = bar_plot, width = 11, height = 8, dpi = 300)

# Bar Plot Facet
# Create a faceted bar plot
facet_bar_plot <- ggplot(prevalence_data, aes(x = MAG, y = Prevalence, fill = Environment)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +  # Replace "Paired" with a dynamic palette
  facet_wrap(~Environment, scales = "free_x", ncol = 1) +  # Facet by Environment
  labs(
    x = "MAG",
    y = "Prevalence (%)",
    title = "Prevalence of Each MAG in Each Environment"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Hide legend to reduce clutter
  )

# Print the faceted bar plot
#print(facet_bar_plot)
output_path <- "./sample_environments/environments/figs/"
ggsave(file.path(output_path, "prevalence_analysis_completed_facet_barplot_only.pdf"), plot = facet_bar_plot, width = 11, height = 8, dpi = 300)


# Heatmap
heatmap_data <- dcast(prevalence_data, MAG ~ Environment, value.var = "Prevalence")
#heatmap_data[is.na(heatmap_data)] <- 0  # Replace NA with 0
heatmap_plot <- ggplot(melt(heatmap_data), aes(x = variable, y = MAG, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(c("darkgrey", "yellow", "red"))(100),
na.value = "white"  # Assign grey80 to NA values
) +
  labs(x = "Environment", y = "MAG", title = "Heatmap of MAG Prevalence Across Environments") +
  theme_classic()
print(heatmap_plot)


# ----------------------------
# 3. Combine All Plots
# ----------------------------

if (!is.null(pie_chart_grid) & !is.null(bar_plot) & !is.null(heatmap_plot)) {
  final_plot <- plot_grid(
    pie_chart_grid,
    facet_bar_plot,
    heatmap_plot,
    labels = c("A", "B", "C"),
    ncol = 1,
    rel_heights = c(0.5, 1.5, 1.5)
  )
  # Save the final plot to a specific path
  output_path <- "./sample_environments/environments/environment_final/figs/"
  ggsave(file.path(output_path, "prevalence_analysis_completed_heatmap_facet_barplot.pdf"), plot = final_plot, width = 12, height = 18, dpi = 300)
} else {
  print("Cannot save the plot, one or more components are missing.")
}


##to add space if necessary
# if (!is.null(pie_chart_grid) & !is.null(facet_bar_plot)) {
#   final_plot <- plot_grid(
#     pie_chart_grid,
#     plot_spacer(),  # Add a spacer plot for space between the two plots
#     facet_bar_plot,
#     labels = c("A", "", "B"),  # Ensure labels are consistent
#     ncol = 1,
#     rel_heights = c(0.5, 0.1, 1.5)  # Adjust heights to control spacing
#   )
#   
#   # Save the final plot to a specific path
#   output_path <- "./sample_environments/environments/environment_final/figs/for_poster_presentation"
#   ggsave(file.path(output_path, "prevalence_analysis_poster_pie_and_barplot.pdf"), plot = final_plot, width = 12, height = 18, dpi = 300)
# } else {
#   print("Cannot save the plot, one or more components are missing.")
# }




###  heatmap plus number of mags
# Load required libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

# Preprocess data for heatmap
heatmap_data <- dcast(prevalence_data, MAG ~ Environment, value.var = "Prevalence")
heatmap_melted <- melt(heatmap_data)

# Heatmap Plot
heatmap_plot <- ggplot(heatmap_melted, aes(x = variable, y = MAG, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colorRampPalette(c("darkgrey", "yellow", "red"))(100),
    na.value = "white",  # Handle missing values
    name = "Prevalence"
  ) +
  labs(x = "Environment", y = "MAG", title = "Heatmap of MAG Prevalence Across Environments") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

# Preprocess data for the small stacked bar plot
stacked_bar_data <- prevalence_data %>%
  group_by(Environment) %>%
  summarise(Unique_MAGs = n_distinct(MAG), .groups = "drop")  # Count unique MAGs per environment

# Simplified stacked bar plot
stacked_bar_plot <- ggplot(stacked_bar_data, aes(x = Environment, y = Unique_MAGs)) +
  geom_bar(stat = "identity", fill = "black") +  # Use a single color for all bars
  labs(
    x = NULL,
    y = NULL,
    title = "Number of Unique MAGs Detected in Each Environment"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"  # Remove the legend
  )

# Combine the heatmap and stacked bar plot
final_plot <- stacked_bar_plot / heatmap_plot + plot_layout(heights = c(1, 4))

# Print the final combined plot
print(final_plot)

#Facet Scatter_plot
# Facet scatter plot by Environment
scatter_plot_facet <- ggplot(filtered_data, aes(x = Breadth, y = Coverage, color = Environment)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_d() +
  facet_wrap(~Environment, scales = "free") +
  labs(
    x = "Breadth",
    y = "Coverage",
    title = "Breadth vs. Coverage Across Environments"
  ) +
  theme_classic()

# Print the faceted plot
print(scatter_plot_facet)


###### make breadth v coverage plots

# Remove rows with missing values for Breadth and Coverage
filtered_data <- data %>%
  filter(!is.na(Breadth) & !is.na(Coverage))

# Create the scatter plot
scatter_plot <- ggplot(filtered_data, aes(x = Breadth, y = Coverage, color = Environment)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_d() +
  labs(
    x = "Breadth",
    y = "Coverage",
    title = "Breadth vs. Coverage Across Environments"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(scatter_plot)

#Facet Scatter_plot
# Facet scatter plot by Environment
scatter_plot_facet <- ggplot(filtered_data, aes(x = Breadth, y = Coverage, color = Environment)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_d() +
  facet_wrap(~Environment, scales = "free") +
  labs(
    x = "Breadth",
    y = "Coverage",
    title = "Breadth vs. Coverage Across Environments"
  ) +
  theme_classic()

# Print the faceted plot
print(scatter_plot_facet)




