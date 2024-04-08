#!/usr/bin/env Rscript

render_qcplots <- function(prefix, nClusters, GMKey){
    #specify the snakemake pipeline working d to knit with
    workd <- getwd()
    output.path <- gsub(".idepth", "_qc.html", normalizePath(paste0(prefix, ".idepth"))) #generate full path of output - brute force because I had issues with relative paths
    
    script.in <- paste0(snakemake@scriptdir, "/qc_dashboard_interactive.Rmd") #get real path of dashboard script
    script.out <- gsub(".Rmd", ".html", paste0(snakemake@scriptdir, "/qc_dashboard_interactive.Rmd")) #get name of future html

    rmarkdown::render(script.in, #knit the markdown file to html
                    params = list(prefix = prefix, nClusters = nClusters, GMKey = GMKey), #pass the path to the QC files that are plotted (via snakemake params)
                    knit_root_dir = workd)  #make sure to knit in the working directory of the snakemake run
    
    #move the default html output to the QC folder. This is an inconvenience of knitr, and 
    #the output.file 
    copy_successful <- file.copy(script.out, output.path)

    # Check if the copy was successful
    if (copy_successful) {
    # If the copy succeeded, delete the original file
    file.remove(script.out)
    } else {
    # If the copy failed, print an error message
    cat("snpArcher: Failed to move the qc dashboard html.\n")
    }
}

render_qcplots(snakemake@params[[1]], snakemake@params[[2]], snakemake@params[[3]])