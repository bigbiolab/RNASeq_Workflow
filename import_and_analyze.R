# Import the data from salmon quant to DESeq2 using tximport
BiocManager::install("tximport")

# Load libraries
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(EnhancedVolcano)


# Get the mapping from transcript IDs to gene symbols 
# What are the columns in the database?
columns(EnsDb.Hsapiens.v86)

# Get the TXID and SYMBOL columns for all entries in database
tx2gene <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                   keys = keys(EnsDb.Hsapiens.v86),
                                   columns = c('TXID', 'SYMBOL'))

# Remove the gene ID column
tx2gene <- dplyr::select(tx2gene, -GENEID)

# Get the quant files and metadata 
# Collect the sample quant files
samples <- list.dirs(
  'Salmon.out/', 
  recursive = FALSE, 
  full.names = FALSE
)

quant_files <- file.path('Salmon.out', samples, 'quant.sf')
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
file.exists(quant_files)  # all should be TRUE

# Set up metadata frame
colData <- data.frame(
  row.names = samples,
  condition = rep(c('untreated', 'dex'), 4)
)

# Compile the tximport counts object and make DESeq dataset
# Get tximport counts object
txi <- tximport(files = quant_files, 
                type = 'salmon',
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = colData,
                                design = ~condition)

# Do DESeq analysis 
# PCA
vsd <- vst(dds)
plotPCA(vsd)

# DEG analysis
dds <- DESeq(dds)

# MA plot 
plotMA(dds)

# get the results 
res <- results(dds)
summary(res)

# Get the results
resdf <- results(dds)
resdf <- as.data.frame(resdf)

