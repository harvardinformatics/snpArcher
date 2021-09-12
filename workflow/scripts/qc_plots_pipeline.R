#!/usr/bin/env Rscript

library(tidyverse)
library(patchwork)
library(plotly)

generate_qcplots <- function(prefix){
  
  # depth plot --------------------------------------------------------------
  
  depth.path <- paste0(prefix, ".idepth")
  df.depth <- read.table(depth.path, header = T)
  
  p.depth <- ggplot(data = df.depth) +
    geom_histogram(aes(x = MEAN_DEPTH), bins = nrow(df.depth)/4,
                   fill = "grey80", color = "black") +
    theme_bw() +
    labs(x = "SNP Depth", y = "Number of individuals")

  
  # PCA plot --------------------------------------------------------------
  pca.path <- paste0(prefix, ".eigenvec")
  tmp.head <- sub("#", "", readLines(pca.path))
  df.pca <- read.table(text = tmp.head, header = TRUE)
  df.val <- read.table(gsub("vec","val", pca.path), header = FALSE)
  df.val$prop <- (df.val$V1 / (sum(df.val$V1))) * 100
  
  #add depth
  df.pca <- left_join(df.pca, df.depth, by = c("IID" = "INDV"))
  
  p.pca <- ggplot(data = df.pca, aes(x = PC1, y = PC2, text = IID, color = MEAN_DEPTH)) + 
    geom_point(alpha = 0.8, shape = 21, fill = "black") + 
    theme_bw() +
    theme(text = element_text(size = 14),
          legend.position = "none") +
    xlab(paste0("PC1", ": ", round(df.val[1,2],1),"% variance")) +
    ylab(paste0("PC2", ": ", round(df.val[2,2],1),"% variance"))
  
  
  # admixture ---------------------------------------------------------------
  
  k2 <- paste0(prefix, ".2.Q")
  k3 <- paste0(prefix, ".3.Q")
  samps <- paste0(prefix, ".fam")
  
  x <- read.table(k2, header = F)
  head(x)
  
  struct_files <- c(k2,k3)
  cat_snmf <- do.call("rbind",lapply(struct_files,
                                 FUN=function(files){
                                   x <- read.table(files, header = F)
                                   names(x) <- gsub("V", "pop", names(x)) #rename ancestral pops
                                   x.samps <- read.table(samps) %>% select(V2) #get sample names from .fam file
                                   x$sampleID <- x.samps$V2 #add sample name to df
                                   x$k <- gsub(".Q","",substr(files, nchar(files)-3+1, nchar(files)))
                                   x.long <- x %>% #pivot longer 
                                     pivot_longer(names_to = "popGroup", values_to = "prob", cols = -c(sampleID, k))
                                   x.long
                                 }))
  
  #nb.cols <- 10
  #myColors <- colorRampPalette(brewer.pal(10, "Set3"))(nb.cols)
  #names(myColors) <- levels(as.factor(cat_snmf$k))
  #colScale <- scale_colour_manual(name = "grp",values = myColors, guide = F)
  #filScale <- scale_fill_manual(name = "grp",values = myColors, guide = F)
  
  #png("output/snmf/admixture_all_k_all_geospiza_big.png", width = 35, height = 25)
  #only ggsave as png works for some odd reason
  
  p.k2 <- cat_snmf %>% 
    filter(k == 2 | k == 3) %>% 
    ggplot(aes(x = factor(sampleID), y = prob, fill = factor(popGroup), text = sampleID)) +
    geom_col(aes(fill = factor(popGroup)), size = 0.1) +
    facet_grid(rows = vars(k), switch = "x", scales = "free", space = "free") +
    theme_minimal() +
    labs(x = "Individuals", y = "Ancestry") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 1)) +
    theme(
      panel.spacing.x = unit(0.0, "lines"),
      axis.text.x = element_text(angle = 90, size = 6),
      panel.grid = element_blank(),
      strip.text.x = element_text(angle = 90, size = 6),
      legend.position = "none"
    )# +
  #scale_fill_manual(name = "grp",values = myColors, guide = F) +xlab(NULL)  + theme(legend.position = "none")

  # SNP quality -------------------------------------------------------------
  
  qual.path <- paste0(prefix, "_snpqc.txt")
  
  df.qual <- read.table(qual.path, na.strings = ".")
  names(df.qual) <- c("CHROM","POS","ID","QUAL","AF","ReadPosRankSum","FS","SOR","MQ","MQRankSum")
  
  # read pos rank sum -------------------------------------------------------

  RPRS.fail <- df.qual %>% filter(ReadPosRankSum < -8) %>% nrow()
  FS.fail <- df.qual %>% filter(FS > 60) %>% nrow()
  SOR.fail <- df.qual %>% filter(SOR > 3) %>% nrow()
  MQ.fail <- df.qual %>% filter(MQ < 40.0) %>% nrow()
  MQRS.fail <- df.qual %>% filter(MQRankSum < -12.5) %>% nrow()
  QUAL.fail <- df.qual %>% filter(QUAL < 30) %>% nrow()
  
  qc.plot <- function(filt.var, cutoff, num.fail){
    x.var <- rlang::sym(filt.var)
    
    if(sum(is.na(df.qual[filt.var])) == nrow(df.qual)){
      p.qc <- ggplot(df.qual)
    }else{
      p.qc <- ggplot(data = df.qual) + 
        geom_freqpoly(aes(x = ! ! x.var), bins = 100) +
        theme_bw() +
        geom_vline(aes(xintercept = cutoff), color = "red") +
        ggtitle(label = paste0("sites filtered = ", num.fail))
    }
    p.qc
  }
  
  p.RPRS <- qc.plot("ReadPosRankSum", -8, RPRS.fail)
  p.FS <- qc.plot("FS", 60, FS.fail)
  p.SOR <- qc.plot("SOR", 3, SOR.fail)
  p.MQ <- qc.plot("MQ",40.0, MQ.fail)
  p.MQRS <- qc.plot("MQRankSum", -12.5, MQRS.fail)
  p.QUAL <- qc.plot("QUAL", 30, QUAL.fail)
  
  layout <- "
  AABB
  AABB
  CCCC
  CCCC
  DEF#
  GH##
  "
  p.pca + p.depth + p.k2 + p.RPRS + p.FS + p.SOR + p.MQ + p.MQRS +
    plot_layout(design = layout)
  ggsave(paste0(prefix, "_qc.pdf"), width = 8.5, height = 11)
  
  # save as interactive html ------------------------------------------------------------

  #make pca interactive
  pl.pca <- ggplotly(p.pca, tooltip=c("text","MEAN_DEPTH", "x","y")) 
  #make admix interactive
  pl.admix <- ggplotly(p.k2, tooltip=c("text", "y"))
  #make first row
  pl.row1 <- subplot(pl.pca, p.depth, titleY = TRUE, titleX = TRUE, margin = 0.06)
  #combine first two plots above admix
  pl.out <- subplot(pl.row1, pl.admix, nrows = 2, margin = 0.08,
                    titleY = TRUE, titleX = TRUE) 
  #save as an html file
  htmlwidgets::saveWidget(as_widget(pl.out), paste0(prefix, "_qc.html"))
  
}

generate_qcplots(snakemake@params[[1]])
