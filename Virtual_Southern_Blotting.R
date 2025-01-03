#!/usr/bin/env Rscript
# ==============================================================================
# Integrated R Script for Enzyme Data
#
# Description:
#   This script combines two functionalities:
#   1) Plotting enzyme data using size-position relationships.
#   2) Converting size columns (e.g., *_H, *_L) to positions and saving results.
#
# Usage (example):
#   Rscript integrated_enzyme_script.R plot \
#       --enzyme-file /path/to/final_distances_with_lengths.csv \
#       --size-position-file /path/to/Souther_coordinate.csv \
#       --plot-output /path/to/enzyme_plot.png
#
#   Rscript integrated_enzyme_script.R convert \
#       --enzyme-file /path/to/final_distances_with_lengths_H.csv \
#       --size-position-file /path/to/Souther_coordinate.csv \
#       --converted-output /path/to/converted_position_data_L.csv
#
# Author:
#   Your Name (your.email@domain.com)
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

# ------------------------------------------------------------------------------
# 1) Plotting Function
# ------------------------------------------------------------------------------
plot_enzyme_data <- function(enzyme_file, size_position_file, plot_output = NULL) {
  # Load enzyme data
  enzyme_data <- read.csv(enzyme_file, stringsAsFactors = FALSE)

  # Load size-position relationship
  size_position <- read.csv(size_position_file, stringsAsFactors = FALSE)

  # Prepare the enzyme data for plotting
  # Assume the columns "IntegrationSite", "Position", and "Chromosome" exist
  # and that all other columns are potential size data
  enzyme_data_long <- pivot_longer(
    data  = enzyme_data,
    cols  = -c("IntegrationSite", "Position", "Chromosome"),
    names_to = "Enzyme",
    values_to = "Size"
  )

  # Function to map size to position based on the given relationship
  map_size_to_position <- function(size) {
    approx(size_position$Size, size_position$Position, xout = size)$y
  }

  # Add a new Position column based on size
  enzyme_data_long$Position <- map_size_to_position(enzyme_data_long$Size)

  # Create the plot
  p <- ggplot(enzyme_data_long, aes(x = Enzyme, y = Position)) +
    geom_point(alpha = 0.5) +
    # Create semi-transparent rectangles for each enzyme
    geom_rect(
      data  = distinct(enzyme_data_long, Enzyme, .keep_all = TRUE),
      aes(xmin = Enzyme, xmax = Enzyme, ymin = -Inf, ymax = Inf),
      alpha = 0.2
    ) +
    scale_y_continuous(
      "Position",
      breaks = size_position$Position,
      labels = size_position$Size
    ) +
    labs(
      title = "Enzyme Sizes",
      x     = "Enzyme",
      y     = "Size"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  # If a plot output is specified, save it; otherwise print to screen
  if (!is.null(plot_output) && nzchar(plot_output)) {
    ggsave(filename = plot_output, plot = p, width = 8, height = 6, dpi = 300)
    message("[INFO] Plot saved to: ", plot_output)
  } else {
    print(p)
    message("[INFO] Plot displayed on screen (no output file specified).")
  }
}

# ------------------------------------------------------------------------------
# 2) Converting Size Columns to Positions
# ------------------------------------------------------------------------------
convert_enzyme_data <- function(enzyme_file, size_position_file, converted_output) {
  # Load enzyme data
  enzyme_data <- read.csv(enzyme_file, stringsAsFactors = FALSE)

  # Load size-position relationship
  size_position <- read.csv(size_position_file, stringsAsFactors = FALSE)

  # Function to map size to position using interpolation
  map_size_to_position <- function(size) {
    approx(size_position$Size, size_position$Position, xout = size)$y
  }

  # Identify columns that match your naming convention (_H or _L) and are numeric
  size_columns <- names(enzyme_data)[
    grepl("(_H|_L)$", names(enzyme_data)) & sapply(enzyme_data, is.numeric)
  ]

  # Convert sizes in those columns to position
  enzyme_data_converted <- enzyme_data
  enzyme_data_converted[size_columns] <- map_df(
    enzyme_data[size_columns],
    map_size_to_position
  )

  # Save updated data
  write.csv(enzyme_data_converted, converted_output, row.names = FALSE)
  message("[INFO] Converted data saved to: ", converted_output)
}

