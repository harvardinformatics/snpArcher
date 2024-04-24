#!/usr/bin/env Rscript

suppressMessages(library(zoo))
suppressMessages(library(patchwork))
suppressMessages(library(gridExtra))
suppressMessages(library(tidyverse))
suppressMessages(library(ggridges))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install("GenomicRanges")
suppressMessages(library(GenomicRanges))

project <- "data/ROH/CCGP/test_qc"

render_rohplots <- function(project){
  #specify the snakemake pipeline working d to knit with
  workd <- getwd()
  output_prefix <- gsub(".rg.roh", "", normalizePath(paste0(project, ".rg.roh"))) #generate full path of output - brute force because I had issues with relative paths
  
  inputFile <- paste0(project, ".rg.roh")
  
  pi1 <- paste0(project, ".1.windowed.pi")
  pi2 <- paste0(project, ".2.windowed.pi")
  
  topF <- paste0(project, "_filtered_top.froh")
  
  # Define a function to load each file into a data frame
  load_file <- function(filename) {
    df <- read.table(filename, header=TRUE)
    df$pi_rm <- zoo::rollmean(df$PI, 10, fill=NA) # apply rollmean to "PI" column
    return(df)
  }
  
  file_list <- list()
  # Use lapply to apply the load_file function to each filename
  file_list <- lapply(paste0(project, ".", 1:10, ".windowed.pi"), load_file)
  
  # Use append to combine the data frames into a single list
  combined_data <- list(file_list)
  
  # get top 2 samples -------------------------------------------------------
  topF <- read.table(topF, col.names = F)
  
  sample1 <- topF[1,1]
  sample2 <- topF[2,1]
  
  # function from DetectRuns ------------------------------------------------
  # not avail on conda so cant put through pipeline
  # https://github.com/bioinformatics-ptp/detectRUNS
  
  subsetBCF <- grep(pattern = "RG", x = readLines(inputFile), 
                    invert = F, value = T)
  
  BCFfinal <- read.table(text = gsub("\t", " ", subsetBCF), 
                         header = F, colClasses = c(rep("character", 3), 
                                                    rep("numeric", 4)), col.names = c("group", "id", 
                                                                                      "chrom", "from", "to", "lengthBps", "nSNP", 
                                                                                      "Quality"))
  BCFfinal$chrom <- gsub("Chr", "", BCFfinal$chrom)
  df_runs = BCFfinal[, c("group", "id", "chrom", "nSNP", 
                         "from", "to", "lengthBps")]
  
  df_runs <- df_runs %>%
    
    # set lengths of runs variable
    mutate(length_bin = ifelse(lengthBps < 100000, "<100kbp",
                               ifelse(lengthBps > 100000 & lengthBps < 500000, "<500kbp",
                                      ifelse(lengthBps > 500000, ">500kbp", NA))))
  
  # -------------------------------------------------------------------------
  # get pi 
  
  #sample 1
  pi1 <- file_list[[1]]
    #read_tsv(pi1, col_names = T)
  roh1 <- df_runs %>% filter(id == sample1) 
  
  pi2 <- file_list[[2]]
    #read_tsv(pi2, col_names = T)
  roh2 <- df_runs %>% filter(id == sample2) 
  
  #rolling mean for pi
  #pi1$pi_rm <- zoo::rollmean(pi1$PI,10,fill=NA)
  #pi2$pi_rm <- zoo::rollmean(pi2$PI,10,fill=NA)
  
  
  # identify largest scaffold -----------------------------------------------
  
  biggest_chr <- pi1 %>% 
    group_by(CHROM) %>% 
    slice_max(order_by = BIN_END, n = 1) %>% 
    arrange(-BIN_END) %>% 
    head(n = 1)
  
  # only chr >1Mbp
  large_chr <- pi1 %>% 
    group_by(CHROM) %>% 
    slice_max(order_by = BIN_END, n = 1) %>% 
    arrange(-BIN_END) %>%
    filter(BIN_END > 1000000) 
  
  # -------------------------------------------------------------------------
  
  # Load necessary packages
  # Set the number of plots per page
  plots_per_page <- 10
  
  # Create a list of individuals to plot
  indiv <-c(topF[1][[1]])
  
  # Create a vector of the chromosome variables to plot
  chromosomes <- large_chr$CHROM
  
  # Create a multipage pdf file to output the plots
  
  pdf(paste0(output_prefix, "_all_indiv_roh.pdf"), width = 8.5, height = 11)
  
  # Loop over the chromosome variables and plot each dataframe
  for (chr in chromosomes) {
    
    # Create an empty list to hold the plots for this chromosome
    plots_list <- list()
    
    # Loop over the data frames and create a plot for each one
    for (i in 1:length(file_list)) {
      
      # Subset the data frame to include only the current chromosome
      df_chr <- file_list[[i]] %>% 
        filter(CHROM == chr)
      
      box_filt <- df_runs %>% filter(id == indiv[i]) %>% 
        filter(chrom == chr)
      
      pi_range <- (max(df_chr$PI) - min(df_chr$PI)) / 50
      
      # Create a plot of the data for this chromosome and add it to the list
      roh_plot <- ggplot(df_chr) + 
        geom_line(aes(x = BIN_START/1000000, y=pi_rm)) +
        geom_rect(data = box_filt, aes(xmin = from/1000000, xmax = to/1000000, 
                                       ymin = -pi_range, ymax =  0, fill = length_bin)) +
        theme_bw() +
        labs(x = NULL, subtitle = (paste0("Sample ", indiv[[i]], ", Chromosome ", chr))) +
        ylab(expression(pi)) +
        scale_fill_manual(values = c("grey80", "blue", "red")) +
        theme(legend.position = "none")
      
      
      
      plots_list[[i]] <- roh_plot
      
      # If we have generated the maximum number of plots per page, create a new page
      if (i %% plots_per_page == 0) {
        grid.arrange(grobs = plots_list, nrow = plots_per_page)
        plots_list <- list()
      }
    }
    
    # If we have any plots left over, create a new page to display them
    if (length(plots_list) > 0) {
      grid.arrange(grobs = plots_list, ncol = plots_per_page)
    }
  }
  
  # Close the pdf file
  dev.off()
  
  
  # build viz ---------------------------------------------------------------
  
  box_filt1 <- roh1 %>% 
    filter(chrom %in% biggest_chr$CHROM)
  
  pi_range <- (max(pi1$PI) - min(pi1$PI)) / 50
  
  title1 <- paste0(basename(project), ": ", sample1, ". Largest scaffold")
  scaff <- paste0(unique(box_filt1$chrom), ": Position (Mbp)")
  
  pA <- pi1 %>% 
    filter(CHROM %in% biggest_chr$CHROM) %>% 
    ggplot() + 
    geom_line(aes(x = BIN_START/1000000, y=pi_rm)) +
    geom_rect(data = box_filt1, aes(xmin = from/1000000, xmax = to/1000000, 
                                    ymin = -pi_range, ymax =  0, fill = length_bin)) +
    theme_bw() +
    labs(x = NULL, title =  title1) +
    ylab(expression(pi)) +
    scale_fill_manual(values = c("grey80", "blue", "red")) +
    theme(legend.position = "bottom")
  
  box_filt2 <- roh2 %>% 
    filter(chrom %in% biggest_chr$CHROM)
  
  pi_range <- (max(pi2$PI) - min(pi2$PI)) / 50
  
  title2 <- paste0(basename(project), ": ", sample2, ". Largest scaffold")
  
  pB <- pi2 %>% 
    filter(CHROM %in% biggest_chr$CHROM) %>% 
    ggplot() + 
    geom_line(aes(x = BIN_START/1000000, y=pi_rm)) +
    geom_rect(data = box_filt2, aes(xmin = from/1000000, xmax = to/1000000, 
                                    ymin = -pi_range, ymax =  0, fill = length_bin)) +
    theme_bw() +
    labs(x = scaff, title =  title2) +
    ylab(expression(pi)) +
    scale_fill_manual(values = c("grey80", "blue", "red")) +
    theme(legend.position = "bottom")
  
  # biggest ROH chunk -------------------------------------------------------
  
  box_largest1 <- roh1 %>% 
    arrange(-lengthBps) %>% 
    head(n = 1)
  
  #scale these plots to the "largest chr" pi
  pi_big <- pi1 %>% 
    filter(CHROM %in% biggest_chr$CHROM) 
  
  title3 <- paste0(basename(project), ": ", sample1, ". Biggest ROH")
  scaff1 <- paste0(unique(box_largest1$chrom), ": Position (Mbp)")
  
  pC <- pi1 %>% 
    filter(CHROM %in% box_largest1$chrom) %>% 
    ggplot() + 
    geom_line(aes(x = BIN_START/1000000, y=pi_rm)) +
    geom_rect(data = box_largest1, aes(xmin = from/1000000, xmax = to/1000000, 
                                       ymin = -pi_range, ymax =  0, fill = length_bin)) +
    theme_bw() +
    labs(x = scaff1, title =  title3) +
    ylab(expression(pi)) +
    scale_fill_manual(values = c( "red")) +
    theme(legend.position = "bottom") +
    ylim(-pi_range, max(pi_big$pi_rm, na.rm = T))
  
  
  box_largest2 <- roh2 %>% 
    arrange(-lengthBps) %>% 
    head(n = 1)
  
  pi_lg2 <- pi2 %>% 
    filter(CHROM %in% box_largest2$chrom)
  
  #scale these plots to the "largest chr" pi
  pi_big2 <- pi2 %>% 
    filter(CHROM %in% biggest_chr$CHROM) 
  
  title4 <- paste0(basename(project), ": ", sample2, ". Biggest ROH")
  scaff2 <- paste0(unique(box_largest2$chrom), ": Position (Mbp)")
  
  pD <- pi2 %>% 
    filter(CHROM %in% box_largest2$chrom) %>% 
    ggplot() + 
    geom_line(aes(x = BIN_START/1000000, y=pi_rm)) +
    geom_rect(data = box_largest2, aes(xmin = from/1000000, xmax = to/1000000, 
                                       ymin = -pi_range, ymax =  0, fill = length_bin)) +
    theme_bw() +
    labs(x = scaff2, title = title4) +
    ylab(expression(pi)) +
    scale_fill_manual(values = c("red")) +
    theme(legend.position = "bottom") +
    ylim(-pi_range, max(pi_big2$pi_rm, na.rm = T))
  
  pA / pB / pC / pD + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  ggsave(paste0(output_prefix, "_pi_roh_top.pdf"), width = 8.5, height = 11)
  
  # check ranges ------------------------------------------------------------
  
  # Create two GRanges objects from the data frames
  gr1 = with(pi1, GRanges(seqnames = CHROM, ranges = IRanges(BIN_START, BIN_END)))
  gr2 = with(roh1, GRanges(seqnames = chrom, ranges = IRanges(from, to)))
  
  # Find overlapping ranges between the two GRanges objects
  overlaps = suppressWarnings(findOverlaps(gr1, gr2))
  
  # Add a column to the first data frame indicating if the range overlaps the ranges in the second data frame
  pi1$overlaps = 0
  pi1$overlaps[queryHits(overlaps)] = 1
  
  pi1$length_bin = NA
  pi1$length_bin[queryHits(overlaps)] = roh1[subjectHits(overlaps), "length_bin"]

  ggplot(pi1, aes(x = PI, y = length_bin)) +
    geom_density_ridges() +
    theme_bw() +
    labs(y = "ROH Length", x = expression(pi), title = sample1)
  ggsave(paste0(output_prefix, "_roh_pi_ridges.pdf"), width = 8, height = 4.5)
  
  # fraction gr 500kb -------------------------------------------------------
  # now calculated earlier in this rule
  
  # approx_length <- pi1 %>% 
  #   group_by(CHROM) %>% 
  #   summarise(max_chr = max(BIN_END)) %>% 
  #   summarise(total_l = sum(max_chr))
  
  # approx_length <- approx_length$total_l
  
  # df_runs %>% 
  #   filter(lengthBps > 500000) %>% 
  #   group_by(id) %>% 
  #   summarise(sum_l = sum(lengthBps)) %>% 
  #   mutate(Froh = sum_l / approx_length) %>% 
  #   write_csv(paste0(output_prefix, "_froh_gr500.csv"))
  
}

render_rohplots(snakemake@params[[1]])