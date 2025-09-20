# Annotating Script
# --- IGNORE ---
# Setup for First Time Only:
# source("setup.R") # Run once to set up the environment

# Setup for All Other Times:
# renv::restore()

# Load Libraries
library(ChIPseeker)
library(ChIPpeakAnno)
library(AnnotationDbi)
library(biomaRt)
library(rtracklayer)

# Databases
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)



library(dplyr)
library(GenomicRanges)

reference_genome <- rtracklayer::import("reference_genome/human/hg19/hg19.fa")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Annotation Example for Peaks
# peaks <- readPeakFile("path_to_peak_file.bed")
# peakAnno <- annotatePeak(peaks, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
# plotAnnoBar(peakAnno)
# plotDistToTSS(peakAnno)

# Annotation Example for Regions
# promoters <- promoters(genes(txdb), upstream=2000, downstream=200)







