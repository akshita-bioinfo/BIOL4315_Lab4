# BIOL4315-Lab4

#loading packages
library(systemPipeRdata)
library(DT)
library(QuasR)
library(Rsubread)

# Ensure required packages are installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("clusterProfiler", "EnhancedVolcano", "biomaRt", "org.At.tair.db", "Rsubread"))

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
# hisat2 -x outputs/hisat2_index/tair10_1_index outputs/processed_fastq/SRR446027_1.fastq.gz_processed.fastq -p 8 -S outputs/sam_files/file1.sam


# creating req directories
dir.create("outputs/bam_sorted_files", recursive = TRUE)
dir.create("outputs/bam_indexed_files", recursive = TRUE)

processed_files <- list.files(
  "outputs/processed_fastq",
  pattern = "_1\\.fastq\\.gz_processed\\.fastq$",
  full.names = TRUE
)
meta_data$FileName1 <- processed_files

# change the metadata paths
#meta_data$FileName1 <- paste0("outputs/processed_fastq/", sub("\\.fastq\\.gz$","", basename(meta_data$FileName1)), ".fastq")
#meta_data$FileName2 <- paste0("outputs/processed_fastq/", basename(meta_data$FileName2), "_processed.fastq")


# ALIGNING THE READS
for(i in 1:nrow(meta_data)){
  
  # retrieve info for sample i(all 18 samples, we are only considering forward reads)
  sample_name <- meta_data$SampleName[i]
  file1 <- meta_data$FileName1[i]
  
# HISAT2
hisat2_log <- tryCatch({
  system2(command = "hisat2", 
          args = c("-x", "outputs/hisat2_index/tair10_1_index/tair10_1_index",
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
index_log <- tryCatch({
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


# READING ALIGNMENT STATS
# get the folder
hisat2_logs_dir <- "outputs/bam_files"

# 1. Get a list of all HISAT2 log files
log_files <- list.files(hisat2_logs_dir, pattern = "_hisat2\\.log$", full.names = TRUE)

if(length(log_files) > 0) {
  percent_aligned <- 1:length(log_files)
  
  for (i in seq_along(percent_aligned)) {
    percent_aligned[i] <- readLines(log_files[i])[length(readLines(log_files[i]))]
  }
  
  align_df <- data.frame(samplename = sort(meta_data$SampleName), percent_aligned)
  align_df <- align_df %>% 
    mutate(percent_aligned = as.numeric(stringr::str_split_i(align_df$percent_aligned, "%", 1))) 
  
  head(align_df)
} else {
  print("No HISAT2 log files found yet. Run the alignment step first.")
}

# Q5: ALIGNMENT STATS
# using DT package to output align_df as a table
datatable(align_df)

# Using ggplot2 to plot a box plot
ggplot(data = align_df, mapping = aes(y = percent_aligned)) + geom_boxplot()


# GENERATING COUNT TABLE
# list bam files
bfiles <- list.files("outputs/bam_sorted_files", pattern = "_sorted.bam$", full.names = TRUE)

# Counting how many reads correspond to each gene
gene_count_list <- Rsubread::featureCounts(
  files = bfiles, 
  annot.ext = "data/GCF_000001735.4_TAIR10.1_genomic.gtf", 
  isGTFAnnotationFile = TRUE, 
  allowMultiOverlap = FALSE,  
  isPairedEnd = FALSE, nthreads = 8,
  minMQS = 10, 
  GTF.featureType = "exon",  
  GTF.attrType = "gene_id" 
)

if(exists("gene_count_list")) {
  glimpse(gene_count_list$counts)[1:5,1:5]
}

if(exists("gene_count_list")) {
  glimpse(gene_count_list$annotation)[1:5,]
}

if(exists("gene_count_list")) {
  glimpse(gene_count_list$stat)[,1:5]
}

# Q6 THE COUNT TABLE
# Assigns $count to a variable named count_table
count_table <- gene_count_list$count
head(count_table)

# Removes the _sorted.bam from the column names
colnames(count_table) <- sub("_sorted\\.bam$", "", colnames(count_table))

# Filter out genes with no reads mapped to them across all samples
count_table <- count_table[rowSums(count_table) > 0, ]
head(count_table)

# DT to output the variable as an interactive table
datatable(count_table)


# DATA PREP FOR DESeq2
if(exists("count_table")) {
  coldata <- meta_data %>% dplyr::select(SampleName,SampleLong,Factor) %>% 
    dplyr::mutate(SampleLong=str_split_i(SampleLong, "\\.",1)) %>% 
    dplyr::rename(condition = SampleLong) %>%
    dplyr::mutate(condition = factor(condition)) %>% 
    dplyr::mutate(Factor = factor(Factor)) 
  
  base::rownames(coldata) <- coldata$SampleName
  coldata <- coldata %>% mutate(SampleName = factor(SampleName))
  coldata$type <- factor(rep("single-read", nrow(coldata)))
  
  coldata <- coldata[base::match(base::colnames(count_table), rownames(coldata)),]
  
  dds1 <- DESeqDataSetFromMatrix(countData = count_table,
                                 colData = coldata,
                                 design = ~ condition)
  
  dds2 <- DESeqDataSetFromMatrix(countData = count_table,
                                 colData = coldata,
                                 design = ~ Factor)
}


# SAMPLE CORRELATION (for dds1)
if(exists("dds1")) {
  d <- cor(assay(rlog(dds1)), method = "spearman")
  hc <- hclust(dist(1 - d))
  
  plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2,
             show.node.label = TRUE, no.margin = TRUE)
}


# ANALYZING DIFFERENTIAL GENE EXPRESSION
if(exists("dds1")) {
  dds1_results <- DESeq(dds1)
  dds2_results <- DESeq(dds2)
  
  res1 <- DESeq2::results(dds1_results)
  res2 <- DESeq2::results(dds2_results)
}


# COMPARING VIR, MOCK, AND AVR BROAD OVERVIEW
if(exists("dds1_results")) {
  res_vir_mock <- DESeq2::results(dds1_results, contrast = c("condition", "Vir", "Mock"), alpha = 0.2)
  res_avr_mock <- DESeq2::results(dds1_results, contrast = c("condition", "Avr", "Mock"), alpha = 0.2)
  res_vir_avr <- DESeq2::results(dds1_results, contrast = c("condition", "Vir", "Avr"), alpha = 0.2)
  
  filter_and_count <- function(res_obj, comparison_name, fc_threshold = 2) {
    res_filtered <- res_obj[!is.na(res_obj$padj) & !is.na(res_obj$log2FoldChange), ]
    sig_genes <- res_filtered[abs(res_filtered$log2FoldChange) >= log2(fc_threshold), ]
    up_regulated <- sum(sig_genes$log2FoldChange > 0)
    down_regulated <- sum(sig_genes$log2FoldChange < 0)
    
    return(data.frame(
      Comparison = comparison_name,
      Up_regulated = up_regulated,
      Down_regulated = down_regulated
    ))
  }
  
  results_summary <- rbind(
    filter_and_count(res_vir_mock, "Vir vs Mock"),
    filter_and_count(res_avr_mock, "Avr vs Mock"),
    filter_and_count(res_vir_avr, "Vir vs Avr")
  )
  
  print(results_summary)
  
  plot_data <- results_summary %>%
    pivot_longer(cols = c(Up_regulated, Down_regulated), 
                 names_to = "Regulation", 
                 values_to = "Count") %>%
    mutate(Regulation = factor(Regulation, levels = c("Up_regulated", "Down_regulated")))
  
  p <- ggplot(plot_data, aes(x = Comparison, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() +  
    labs(
      title = "Differentially Expressed Genes by Comparison",
      subtitle = "Fold Change >= 2, alpha = 0.2",
      x = "Comparison",
      y = "Number of Genes",
      fill = "Regulation"
    ) +
    theme_minimal() 
  print(p)
}


# QUESTION 7
# get the comp pairs
comp <- systemPipeR::readComp(system.file("extdata/param/targetsPE.txt", package="systemPipeRdata"))
comp[[1]]

# DGE analysis
if(exists("dds2_results")) {
  res_M1_A1 <- DESeq2::results(dds2_results, contrast = c("Factor", "M1", "A1"), alpha = 0.2)
  res_M1_V1 <- DESeq2::results(dds2_results, contrast = c("Factor", "M1", "V1"), alpha = 0.2)
  res_A1_V1 <- DESeq2::results(dds2_results, contrast = c("Factor", "A1", "V1"), alpha = 0.2)
  res_M6_A6 <- DESeq2::results(dds2_results, contrast = c("Factor", "M6", "A6"), alpha = 0.2)
  res_M6_V6 <- DESeq2::results(dds2_results, contrast = c("Factor", "M6", "V6"), alpha = 0.2)
  res_A6_V6 <- DESeq2::results(dds2_results, contrast = c("Factor", "A6", "V6"), alpha = 0.2)
  res_M12_A12 <- DESeq2::results(dds2_results, contrast = c("Factor", "M12", "A12"), alpha = 0.2)
  res_M12_V12 <- DESeq2::results(dds2_results, contrast = c("Factor", "M12", "V12"), alpha = 0.2)
  res_A12_V12 <- DESeq2::results(dds2_results, contrast = c("Factor", "A12", "V12"), alpha = 0.2)
  
  filter_and_count <- function(res_obj, comparison_name, fc_threshold = 2) {
    res_filtered <- res_obj[!is.na(res_obj$padj) & !is.na(res_obj$log2FoldChange), ]
    sig_genes <- res_filtered[abs(res_filtered$log2FoldChange) >= log2(fc_threshold), ]
    up_regulated <- sum(sig_genes$log2FoldChange > 0)
    down_regulated <- sum(sig_genes$log2FoldChange < 0)
    
    return(data.frame(
      Comparison = comparison_name,
      Up_regulated = up_regulated,
      Down_regulated = down_regulated
    ))
  }
  
  results_summary <- rbind(
    filter_and_count(res_M1_A1, "M1 VS A1"),
    filter_and_count(res_M1_V1, "M1 VS V1"),
    filter_and_count(res_A1_V1, "A1 VS V1"),
    filter_and_count(res_M6_A6, "M6 VS A6"),
    filter_and_count(res_M6_V6, "M6 VS V6"),
    filter_and_count(res_A6_V6, "A6 VS V6"),
    filter_and_count(res_M12_A12, "M12 VS A12"),
    filter_and_count(res_M12_V12, "M12 VS V12"),
    filter_and_count(res_A12_V12, "A12 VS V12")
   
  )
  
  print(results_summary)
  
  plot_data <- results_summary %>%
    pivot_longer(cols = c(Up_regulated, Down_regulated), 
                 names_to = "Regulation", 
                 values_to = "Count") %>%
    mutate(Regulation = factor(Regulation, levels = c("Up_regulated", "Down_regulated")))
  
  p <- ggplot(plot_data, aes(x = Comparison, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() +  
    labs(
      title = "Differentially Expressed Genes by Comparison",
      subtitle = "Fold Change >= 2, alpha = 0.2",
      x = "Comparison",
      y = "Number of Genes",
      fill = "Regulation"
    ) +
    theme_minimal() 
  print(p)
}

# QUESTION 8: COMPARE HEATMAP VS SPEARMAN CLUSTERING

if(exists("dds2")) {
  d <- cor(assay(rlog(dds2)), method = "spearman")
  hc <- hclust(dist(1 - d))
  
  plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2,
             show.node.label = TRUE, no.margin = TRUE)
}

# CODE FOR DESCRIPTIVE VARIBALE
# Connect to Ensembl Plants BioMart and download TAIR gene descriptions  
m <- biomaRt::useMart("plants_mart",                                     
                      dataset = "athaliana_eg_gene",                     
                      host = "https://plants.ensembl.org")               

desc <- AnnotationDbi::select(org.At.tair.db,                            
                              keys = rownames(res_vir_mock),             
                              columns = c("GENENAME"),                   
                              keytype = "TAIR") %>%                      
  dplyr::rename(gene_id = TAIR, description = GENENAME) %>%              
  dplyr::distinct(gene_id, .keep_all = TRUE)

# desc <- biomaRt::getBM(attributes = c("tair_locus", "description"), mart = m)                                                                       
# desc <- desc[!duplicated(desc[, 1]), ]                                   
# desc <- desc %>% dplyr::rename(gene_id = tair_locus)

# ADDING GENE DESCRIPTIONS AND GETTING SPECIFIC WITH VOLACANO PLOTS
if(exists("res_vir_mock") && exists("desc")) {
  annotate_results <- function(res_obj, desc_df) {
    res_df <- as.data.frame(res_obj)
    res_df$gene_id <- rownames(res_df)
    res_df <- left_join(res_df, desc_df, by = "gene_id") %>%
      mutate(description = str_split_i(description,"\\[",1)) 
    return(res_df)
  }
  
  res_vir_mock_annot <- annotate_results(res_vir_mock, desc)
  res_avr_mock_annot <- annotate_results(res_avr_mock, desc)
  res_vir_avr_annot <- annotate_results(res_vir_avr, desc)
  
  volcano1 <- EnhancedVolcano(res_vir_mock_annot,
                              lab = res_vir_mock_annot$description,
                              x = 'log2FoldChange',
                              y = 'pvalue',
                              title = 'Vir vs Mock',
                              pCutoff = 0.05,           
                              FCcutoff = 1.0,
                              pointSize = 4.0,
                              labSize = 4.0,
                              drawConnectors = TRUE)
  print(volcano1)
}


# GENE ONTOLOGY (GO) Enrichment
# Select significant UP-regulated genes from Avirulent vs Mock
# Because we are working with a downsampled toy dataset, we use a relaxed pvalue cutoff
sig_avr_up <- res_avr_mock %>%
  as.data.frame() %>%
  dplyr::filter(pvalue < 0.05 & log2FoldChange > 1) %>%
  rownames()

# Run the enrichment using TAIR IDs
ego_avr <- enrichGO(gene          = sig_avr_up,
                    OrgDb         = org.At.tair.db,
                    keyType       = "TAIR",
                    ont           = "BP", # Biological Process
                    pAdjustMethod = "none",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2)

# Visualize with a dotplot
enrichplot::dotplot(ego_avr, showCategory=20) + 
  ggplot2::ggtitle("GO Enrichment: Avirulent Response (Up-regulated)")