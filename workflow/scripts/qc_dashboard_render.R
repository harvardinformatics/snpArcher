#!/usr/bin/env Rscript

render_qcplots <- function(prefix){
    #specify the snakemake pipeline working d to knit with
    workd <- getwd()
    output.path <- gsub(".idepth", "_qc.html", normalizePath(paste0(prefix, ".idepth"))) #generate full path of output - brute force because I had issues with relative paths

    rmarkdown::render('../workflow/scripts/qc_dashboard_interactive.Rmd', #knit the markdown file to html
                    params = list(prefix = prefix), #pass the path to the QC files that are plotted (via snakemake params)
                    knit_root_dir = workd)  #make sure to knit in the working directory of the snakemake run
    
    #move the default html output to the QC folder. This is an inconvenience of knitr, and 
    #the output.file 
    file.rename(from = "../workflow/scripts/qc_dashboard_interactive.html",
          to = output.path)
}

render_qcplots(snakemake@params[[1]])
