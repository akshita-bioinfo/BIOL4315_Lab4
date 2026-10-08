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


# ACCESSING RNAseq DATA AND SAMPLE DATA                    
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
      

# READ PROCESSING
# consider pairs for all 36 files
all_pairs <- rep(1:18, each = 2)

# rqc run
qc_results <- rqc(
  path = fq_path, 
  pattern = ".fastq.gz", 
  pair = all_pairs, 
  openBrowser = FALSE,
  outdir = "outputs")

# Generate 3 per-cycle Q-score box plots for files 1-12, 13-24, 25-36.
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


# QuasR TRIMMING
# creating folder
dir.create("outputs/processed_fastq", recursive = TRUE)

# iterate over every sample in meta_data table to complete QuasR trimming
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


# ALIGNMENTS

#create the hisat2_index directory
dir.create("outputs/hisat2_index", recursive = TRUE)

# variable to store ref genome
at_genome <- "data/GCF_000001735.4_TAIR10.1_genomic.fna"

# create tair10_1_index directory
dir.create("outputs/hisat2_index/tair10_1_index", recursive = TRUE)

#Use system2 to run hisat2 (terminal cmd) from within R
tryCatch({
  system2(command = "hisat2-build", 
          args = c("-p","8", at_genome, "outputs/hisat2_index/tair10_1_index"),
          stdout = TRUE, stderr = TRUE)
}, error = function(e) {
  paste("hisat2-build", "indexing failed with error:", e$message)
})


# creating outputs folders
dir.create("outputs/sam_files", recursive = TRUE)
dir.create("outputs/bam_files", recursive = TRUE)

# gunzip all processed fastq.gz files
# Gather all 36 processed .gz files
# processed_files <- list.files(path = "outputs/processed_fastq", pattern = "\\.fastq\\.gz_processed\\.fastq\\gz$", full.names = TRUE)
# update metadata table
# meta_data$FileName1 <- list.files(path = "outputs/processed_fastq", pattern = "_1\\.fastq\\.gz_processed\\.fastq$", full.names = TRUE)
# meta_data$FileName2 <- list.files(path = "outputs/processed_fastq", pattern = "_2\\.fastq\\.gz_processed\\.fastq$", full.names = TRUE)
# hisat2 -x outputs/hisat2_index/tair10_1_index outputs/processed_fastq/SRR446027_1.fastq.gz_processed.fastq -p 8 -S outputs/sam_files/file1.sam


# creating req directories
dir.create("outputs/bam_sorted_files", recursive = TRUE)
dir.create("outputs/bam_indexed_files", recursive = TRUE)

# ALIGNING THE READS
for(i in 1:nrow(meta_data)){
  
  # retrieve info for sample i(all 18 samples, we are only considering forward reads)
  sample_name <- meta_data$SampleName[i]
  file1 <- meta_data$FileName1[i]
  
# HISAT2
hisat2_log <- tryCatch({
  system2(command = "hisat2", 
          args = c("-x", "outputs/hisat2_index/tair10_1_index",
                   "-U", file1,
                   "-p", "8" ,
                   "-S",
                   paste0("outputs/sam_files/", sample_name, ".sam")),
          stdout = TRUE, 
          stderr = TRUE)
}, error = function (e) {
  paste("hisat2", "alignment failed with error:", e$message)
})

# WRITE HISAT2 LOG
writeLines(
  hisat2_log, 
  file.path("outputs/bam_files",
  paste0(sample_name, "_hisat2.log")
  )
)

# SAMTOOLS: VIEW BAM FILES
view_log <- tryCatch({
  system2(command = "samtools",
          args = c("view", 
                   paste0("outputs/sam_files/", sample_name, ".sam"),
                   "-b", # output bam format
                   "-o",
                   paste0("outputs/bam_files/", sample_name, ".bam")),
          stdout = TRUE, 
          stderr = TRUE)
}, error = function (e) {
  paste("samtools", "bam coversion failed with error:", e$message)
})

# SAMTOOLS: SORT BAM FILES
sort_log <- tryCatch({
  system2(command = "samtools",
          args = c("sort", 
                   paste0("outputs/bam_files/", sample_name, ".bam"),
                  "-o",
                  paste0("outputs/bam_sorted_files/", sample_name, "_sorted.bam")),
                  stdout = TRUE, 
                  stderr = TRUE)
}, error = function (e) {
  paste("samtools", "sorting failed with error:", e$message)
})

# SAMTOOLS: INDEXING SORTED BAM FILES
index-log <- tryCatch({
  system2(command="samtools",
          args = c("index", 
                    paste0("outputs/bam_sorted_files/", sample_name, "_sorted.bam"),
                    "-o",
                    paste0("outputs/bam_indexed_files/", sample_name, ".bai")),
                    stdout = TRUE, 
                    stderr = TRUE)
}, error = function (e) {
  paste("samtools", "indexing failed with error:", e$message)
})
 
}




