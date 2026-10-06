# BIOL4315-Lab4

#loading packages
library(systemPipeRdata)
library(DT)
library(QuasR)

# Ensure required packages are installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("clusterProfiler", "EnhancedVolcano", "biomaRt", "org.At.tair.db"))

lapply(c(
  "docopt", "DT", "pheatmap", "GenomicFeatures", "DESeq2",
  "edgeR", "systemPipeR", "systemPipeRdata", "BiocStyle", "GO.db", "dplyr",
  "tidyr", "stringr", "Rqc", "QuasR", "ape", "clusterProfiler", "biomaRt", "EnhancedVolcano", "org.At.tair.db"
), require, character.only = TRUE)

# e.g.1 get the path 
fq_path <- systemPipeRdata::pathList()$fastqdir

# e.g.2 list the files 
fq_files <- list.files(fq_path)
head(fq_files)


# read the table
meta_data <- read.table(
  system.file("extdata/param/targetsPE.txt", package="systemPipeRdata"), 
  header = TRUE)
head(meta_data)

# checking for path from package
# actual_path <- systemPipeRdata::pathList()$fastqdir
# actual_path <- paste0(system.file("extdata", package="systemPipeRdata"), "/")
# print(actual_path)



# code chunk to replace ./data/ with actual path
meta_data$FileName1 <- gsub("./data/", fq_path, meta_data$FileName1, fixed = TRUE)
meta_data$FileName2 <- gsub("./data/", fq_path, meta_data$FileName2, fixed = TRUE)

# verify
head(meta_data[, c("FileName1", "FileName2")])

# output interactive table
datatable(meta_data, 
          options = list(
            pagelength = 6, 
            scrollX = TRUE, 
            autoWidth = TRUE
          ),
          caption = 'Table1: Updated Metadata Table')
      

# Generate 3 per-cycle Q-score box plots for files 1-12, 13-24, 25-36.
# consider pairs for all 36 files
all_pairs <- rep(1:18, each = 2)

# rqc run
qc_results <- rqc(
  path = fq_path, 
  pattern = ".fastq.gz", 
  pair = all_pairs, 
  openBrowser = FALSE)

# Q-score box plots
# plot1
rqcCycleQualityBoxPlot(qc_results[1:12])

# plot2
rqcCycleQualityBoxPlot(qc_results[13:24])

# plot3
rqcCycleQualityBoxPlot(qc_results[25:36])

# base call frequency plots
# plot1
rqcCycleBaseCallsLinePlot(qc_results[1:12])

# plot2
rqcCycleBaseCallsLinePlot(qc_results[13:24])

# plot3
rqcCycleBaseCallsLinePlot(qc_results[25:36])


dir.create("outputs/processed_fastq", recursive = TRUE)

for (i in seq_along(meta_data$FileName1)) {
  
  fastqfiles <- c(meta_data$FileName1[i], meta_data$FileName2[i])
  l <- length(stringr::str_split_1(meta_data$FileName1[i], "\\/"))
  
  outfiles <- c(paste0("outputs/processed_fastq/", stringr::str_split_1(meta_data$FileName1[i], "\\/")[l], "_processed.fastq.gz"),
                paste0("outputs/processed_fastq/", stringr::str_split_1(meta_data$FileName2[i], "\\/")[l], "_processed.fastq.gz"))
  
  QuasR::preprocessReads(fastqfiles, outfiles,
                         nBases=1,
                         truncateEndBases=3,
                         Lpattern="GCCCGGGTAA",
                         minLength=40)
  
  meta_data$FileName1[i] <- outfiles[1]
  meta_data$FileName2[i] <- outfiles[2]
}
