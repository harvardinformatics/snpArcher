#!/usr/bin/env Rscript

## Script to clean up messy sample data for snpArcher QC

# Check if the tidyverse package is installed
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  # Install the tidyverse package if it is not installed
  install.packages("tidyverse")
}

# Check if the parzer package is installed
if (!requireNamespace("parzer", quietly = TRUE)) {
  # Install the tidyverse package if it is not installed
  install.packages("parzer")
}

# Load the tidyverse and parzer package
library(tidyverse, quietly=TRUE, warn.conflicts = FALSE)
library(parzer, quietly=TRUE, warn.conflicts = FALSE)

add_lat_long <- function(tbl) {
  # Check if the tibble has a "Latitude" and "Longitude" columns
  if ("Latitude" %in% colnames(tbl) && "Longitude" %in% colnames(tbl)) {
    # If so, create a new "lat" and "long" columns using the values from the "Latitude" and "Longitude" columns
    tbl <- tbl %>%
      mutate(lat = parse_lat(Latitude),
             long = parse_lon(Longitude))
  } else if ("lat_lon" %in% colnames(tbl)) {
    # If the tibble has a "lat_lon" column, create a new "lat" and "long" columns using the values from the "lat_lon" column
    tbl <- tbl %>%
      mutate(lat_lon_clean = str_replace_all(tbl$lat_lon, "[\r\n\"\'\ ́]", " ") %>% 
               str_trim() %>%
               str_replace_all("(.)\\s+([:alpha:])", "\1 , \2")) %>%
      mutate(lat = parse_llstr(lat_lon_clean)$lat,
             long =  parse_llstr(lat_lon_clean)$lon)
  } else {
    tbl <- tbl %>%
      mutate(lat = NA,
             long = NA)
  }
  # Return the updated tibble
  return(tbl)
}

remove_columns_with_only_na <- function(tbl) {
  # Get the names of all columns in the tibble
  cols <- colnames(tbl)
  
  # Filter out columns that contain only NA values
  cols_to_keep <- cols[!sapply(tbl[, cols], function(x) all(is.na(x)))]
  
  # Return the original tibble with only the columns that do not contain only NA values
  tbl[, cols_to_keep, drop = FALSE]
}

# Check that the correct number of command line arguments were provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  stop("Usage: remove_columns_with_only_na.R <input_file> <output_file>")
}

# Get the input and output file names from the command line arguments
input_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

# Read the input file as a tibble
tbl <- read_csv(input_file, col_names=TRUE, show_col_types = FALSE)

# Process and output
tbl %>% 
  remove_columns_with_only_na() %>% 
  add_lat_long() %>% 
  select(BioSample,lat,long,refGenome) %>% 
  arrange(BioSample) %>% 
  distinct(BioSample,lat,long,.keep_all=TRUE) %>%
  remove_columns_with_only_na() %>% write_csv(output_file)
