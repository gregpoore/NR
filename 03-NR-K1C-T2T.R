# 03-NR-K1C-T2T.R
# Author: Greg Poore
# Date: Feb 20, 2024
# Purpose: Analyze TCGA human reads input into Kraken

#-------------------------------#
# Load dependencies
require(doMC)
require(plyr)
require(dplyr)
require(tibble)
require(reshape2)
require(phyloseq)
require(biomformat)

numCores <- detectCores()
registerDoMC(cores=numCores)

# Create figure folder
plotPrefix <- "K1C-T2T"
if(!( dir.exists( file.path( paste0("Figures/",plotPrefix) )))){
  dir.create(file.path( paste0("Figures/",plotPrefix) ))
}
if(!( dir.exists( file.path( paste0("Interim_data/",plotPrefix) )))){
  dir.create(file.path( paste0("Interim_data/",plotPrefix) ))
}

#----------------------------------------------------------#
# Import TCGA data
#----------------------------------------------------------#
# Load original data
load("Input_data/tcgaVbDataAndMetadataAndSNM.RData", verbose = TRUE)

# Load count data
countsVb_BIOM <- read_biom(biom_file = "Input_data/against-tcga-kraken1-processed/t2t-only-pese/genus_dbtcga_without_580_and_29_bac_and_vir_NONE_merged.biom")
countsVb <- data.frame(t(as(biom_data(countsVb_BIOM),"matrix")))

countsVbSampleIDs <- rownames(countsVb)
countsVbSampleIDs2 <- gsub("\\.filtered\\.R1.+","",countsVbSampleIDs)
countsVbSampleIDs3 <- gsub("\\.R1\\.trimmed.+","",countsVbSampleIDs2)
countsVbSampleIDs4 <- gsub("\\.R1\\.trimmed$","",countsVbSampleIDs3)
rownames(countsVb) <- countsVbSampleIDs4

# Load metadata
load("Input_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData", verbose = TRUE)
# metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts
metaMyco <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  select(-bam_mapped_reads, -bam_unmapped_reads, -bam_ratio_unmapped)
metaData <- metaMyco %>% rownames_to_column("mycoIDs")
metaDataFilenames <- gsub("\\.filtered\\.$","",metaData$run_prefix)
rownames(metaData) <- metaDataFilenames

# Check overlap both ways
sum(metaDataFilenames %in% rownames(countsVb)) # 15512
sum(rownames(countsVb) %in% metaDataFilenames) # 15512

# Reorder metadata
metaDataOrd <- metaData[rownames(countsVb),]

# Use myco sample IDs for conciseness
metaDataOrdFor <- metaDataOrd %>%
  rownames_to_column("krIDs") %>%
  column_to_rownames("mycoIDs")

countsVb2 <- countsVb
identical(rownames(countsVb2), metaDataOrdFor$krIDs) # TRUE
rownames(countsVb2) <- rownames(metaDataOrdFor)

# Check final results
head(metaDataOrdFor)
countsVb2[1:3,1:3]
identical(rownames(metaDataOrdFor),rownames(countsVb2)) # TRUE

# Save as clean tables
metaDataFinal <- metaDataOrdFor
countsVbFinal <- countsVb2

# Check for any zero-sum samples and remove -- none are removed
metaDataFinalZeroSum <- metaDataFinal[which(rowSums(countsVbFinal)==0),]
table(metaDataFinalZeroSum$experimental_strategy) # 0

metaDataFinalNonzero <- droplevels(metaDataFinal[which(rowSums(countsVbFinal)!=0),])
countsVbFinalNonzero <- countsVbFinal[which(rowSums(countsVbFinal)!=0),]
dim(metaDataFinalNonzero) # 15512    41
dim(countsVbFinal) # 15512  1874

#----------------------------------------------------------#
# Import WIS data summarized at genus level
#----------------------------------------------------------#
load("Input_data/wis-bacteria-fungi-genera-species-bio-24July22.RData", verbose = TRUE)

wisGenera <- data.frame(tax_table(psWzBacteriaAndFungi_genus_Bio))
wisGeneraBact <- wisGenera %>% filter(kingdom == "Bacteria") %>% droplevels()
wisGeneraBactKnown <- wisGeneraBact %>% filter(!grepl("Unknown",genus)) %>% droplevels()
length(unique(wisGeneraBactKnown$genus)) # 216

wisGeneraBactKnownUnique <- unique(wisGeneraBactKnown$genus)

#----------------------------------------------------------#
# Intersect Kraken-derived TCGA data and WIS genera
#----------------------------------------------------------#

krakenGenera <- gsub("k__.+\\.g__","",colnames(countsVbFinalNonzero))

# Testing special cases
grep("Shigella",krakenGenera, value = TRUE)
grep("Escherichia",krakenGenera, value = TRUE)
grep("Clostridium",krakenGenera, value = TRUE)

# Notes on merging taxa @ genus level:
# - WIS combined "Escherichia/Shigella" and both are found in krakenGenera --> allow both
# - WIS has several versions of Clostridium but Kraken only has one --> allow one

wisGeneraBactKnownUniqueRev <- c(wisGeneraBactKnownUnique,
                                 "Escherichia",
                                 "Shigella")

# Check overlap number
sum(krakenGenera %in% wisGeneraBactKnownUniqueRev) # 184
sum(wisGeneraBactKnownUniqueRev %in% krakenGenera) # 184

# Intersect features
krakenGeneraWISInteresected <- intersect(krakenGenera, wisGeneraBactKnownUniqueRev)

# Create new data frames of intersected dataframes
tcgaGenusKrakenAll <- countsVbFinalNonzero
colnames(tcgaGenusKrakenAll) <- gsub("k__.+\\.g__","",colnames(tcgaGenusKrakenAll))

countsVbFinalNonzeroWIS <- tcgaGenusKrakenAll[, krakenGeneraWISInteresected]

dim(countsVbFinalNonzeroWIS) # 15512   184

# - No samples have zero counts after the feature subsetting
sum(rowSums(countsVbFinalNonzeroWIS)==0) # 0

# save(countsVbFinalNonzeroWIS,
#      metaDataFinalNonzero,
#      file = "Interim_data/tcga-wis-overlapping-data-and-metadata-subset-24July22.RData")

#----------------------------------------------------------#
# Intersect Kraken-derived TCGA data and metagenome Bin genera
#----------------------------------------------------------#

binTumorTaxa <- read.csv("Input_data/bins-method-five-tumor-only/bin-taxonomy.tsv",
                         sep = "\t", stringsAsFactors = FALSE) %>%
  rename("binID" = "X.input_bin",
         "phylophlanString" = "X.u.k._.S.G.F.GBid.taxa_level.taxonomy.avg_dist") %>%
  tidyr::separate("phylophlanString", c("sgbID","sgbHitTaxaLevel","estimated_taxonomy","avgDist"),
                  sep = "\\:", remove = TRUE) %>%
  mutate(estimated_taxonomy = gsub("k__|p__|c__|o__|f__|g__|s__|t__","",estimated_taxonomy)) %>%
  tidyr::separate("estimated_taxonomy", c("domain","phylum","class","order","family","genus","species","strain"),
                  sep = "\\|", remove = TRUE)

binTumorStats <- read.csv("Input_data/bins-method-five-tumor-only/bin-stats.tsv",
                          sep = "\t", stringsAsFactors = FALSE) %>%
  rename("binID" = "X")

sum(binTumorTaxa$binID %in% binTumorStats$binID) # 759
sum(binTumorStats$binID %in% binTumorTaxa$binID) # 759

binTumorTaxaAndStats <- binTumorTaxa %>%
  left_join(binTumorStats, by = "binID")

binTumorTaxaAndStatsKnown <- binTumorTaxaAndStats %>%
  filter(!grepl("^GGB.+",species))

binTumorGeneraKnownUnique <- sort(unique(binTumorTaxaAndStatsKnown$genus))
length(binTumorGeneraKnownUnique) # 179

# Create new data frames of intersected dataframes
tcgaGenusKrakenAll <- countsVbFinalNonzero
colnames(tcgaGenusKrakenAll) <- gsub("k__.+\\.g__","",colnames(tcgaGenusKrakenAll))

countsVbFinalNonzeroBins <- tcgaGenusKrakenAll[, intersect(colnames(tcgaGenusKrakenAll),
                                                           binTumorGeneraKnownUnique)]
dim(countsVbFinalNonzeroBins) # 15512   132

# - No samples have zero counts after the feature subsetting
sum(rowSums(countsVbFinalNonzeroBins)==0) # 0

# Find overlap of WIS and Bins
binWISOverlap <- intersect(wisGeneraBactKnownUniqueRev,
                           binTumorGeneraKnownUnique)
length(binWISOverlap) # 69 (31.7% of 218)

countsVbFinalNonzeroBinsWIS <- tcgaGenusKrakenAll[, intersect(krakenGeneraWISInteresected,
                                                              binTumorGeneraKnownUnique)]
dim(countsVbFinalNonzeroBinsWIS) # 15512   69

#----------------------------------------------------------#
# Subset to QC samples for batch correction
#----------------------------------------------------------#

table(metaDataFinalNonzero$data_submitting_center_label)
table(metaDataFinalNonzero$sample_type)
table(metaDataFinalNonzero$cgc_platform)
table(metaDataFinalNonzero$experimental_strategy)

# Intersect samples with original QC'd metadata
metaDataFinalNonzeroQC <- metaDataFinalNonzero %>%
  filter(knightlabID %in% rownames(metadataSamplesAllQC)) %>%
  droplevels()
dim(metaDataFinalNonzeroQC) # 15180    41

# Subset count data
countsVbFinalNonzeroQC <- countsVbFinalNonzero[rownames(metaDataFinalNonzeroQC),]
countsVbFinalNonzeroQCWIS <- countsVbFinalNonzeroWIS[rownames(metaDataFinalNonzeroQC),]
countsVbFinalNonzeroQCBins <- countsVbFinalNonzeroBins[rownames(metaDataFinalNonzeroQC),]
countsVbFinalNonzeroQCBinsWIS <- countsVbFinalNonzeroBinsWIS[rownames(metaDataFinalNonzeroQC),]

save(countsVbFinalNonzeroQC,
     countsVbFinalNonzeroQCWIS,
     countsVbFinalNonzeroQCBins,
     metaDataFinalNonzeroQC,
     file = paste0("Interim_data/",plotPrefix,"/data_raw_tcga_full_wis_bins_features_subset_23Feb24.RData"))

#----------------------------------------------------#
# Perform batch correction 
#----------------------------------------------------#

## VSNM -- Full data
source("00-functions.R")
# Note: covDesignNorm had dimensions: 17625   211
batchCorrectedVSNMData <- vsnmFunctionTCGA(qcData = countsVbFinalNonzeroQC,
                                           qcMetadata = metaDataFinalNonzeroQC)
voomDataGenusKrakenQCFilt <- batchCorrectedVSNMData$vdge_dataE
vsnmDataGenusKrakenQCFilt <- batchCorrectedVSNMData$snmData

## VSNM -- WIS feature subset data
source("00-functions.R")
# Note: covDesignNorm had dimensions: 15180    18
batchCorrectedVSNMDataWIS <- vsnmFunctionTCGA(qcData = countsVbFinalNonzeroQCWIS,
                                              qcMetadata = metaDataFinalNonzeroQC)
voomDataGenusKrakenQCFiltWIS <- batchCorrectedVSNMDataWIS$vdge_dataE
vsnmDataGenusKrakenQCFiltWIS <- batchCorrectedVSNMDataWIS$snmData

## VSNM -- Bins feature subset data
source("00-functions.R")
# Note: covDesignNorm had dimensions: 15180    18
batchCorrectedVSNMDataBins <- vsnmFunctionTCGA(qcData = countsVbFinalNonzeroQCBins,
                                               qcMetadata = metaDataFinalNonzeroQC)
voomDataGenusKrakenQCFiltBins <- batchCorrectedVSNMDataBins$vdge_dataE
vsnmDataGenusKrakenQCFiltBins <- batchCorrectedVSNMDataBins$snmData

# # Not currently converging
# ## VSNM -- BinsWIS feature subset data
# source("00-functions.R")
# # Note: covDesignNorm had dimensions: 15180    18
# batchCorrectedVSNMDataBinsWIS <- vsnmFunctionTCGA(qcData = countsVbFinalNonzeroQCBinsWIS,
#                                                qcMetadata = metaDataFinalNonzeroQC)
# voomDataGenusKrakenQCFiltBinsWIS <- batchCorrectedVSNMDataBinsWIS$vdge_dataE
# vsnmDataGenusKrakenQCFiltBinsWIS <- batchCorrectedVSNMDataBinsWIS$snmData

save(countsVbFinalNonzeroQC,
     voomDataGenusKrakenQCFilt,
     vsnmDataGenusKrakenQCFilt,
     countsVbFinalNonzeroQCWIS,
     voomDataGenusKrakenQCFiltWIS,
     vsnmDataGenusKrakenQCFiltWIS,
     countsVbFinalNonzeroQCBins,
     voomDataGenusKrakenQCFiltBins,
     vsnmDataGenusKrakenQCFiltBins,
     # countsVbFinalNonzeroQCBinsWIS,
     # voomDataGenusKrakenQCFiltBinsWIS,
     # vsnmDataGenusKrakenQCFiltBinsWIS,
     metaDataFinalNonzeroQC,
     file = paste0("Interim_data/",plotPrefix,"/data_vsnm_tcga_full_wis_bins_features_subset_20Feb24.RData"))

#----------------------------------------------------#
# Subset to per-center data
#----------------------------------------------------#

metaDataFinalNonzeroQC %>% count(cgc_platform) # HiSeq accounts for 98.2% of samples (14913/15179)
metaDataFinalNonzeroQC %>% count(cgc_platform, data_submitting_center_label)

#-----------------Subset to Illumina HiSeq-----------------#
metaDataFinalNonzeroQC_HiSeq <- metaDataFinalNonzeroQC %>%
  filter(cgc_platform == "Illumina HiSeq") %>% droplevels()

countsVbFinalNonzeroQC_HiSeq <- countsVbFinalNonzeroQC[rownames(metaDataFinalNonzeroQC_HiSeq),]
countsVbFinalNonzeroQCWIS_HiSeq <- countsVbFinalNonzeroQCWIS[rownames(metaDataFinalNonzeroQC_HiSeq),]
countsVbFinalNonzeroQCBins_HiSeq <- countsVbFinalNonzeroQCBins[rownames(metaDataFinalNonzeroQC_HiSeq),]
countsVbFinalNonzeroQCBinsWIS_HiSeq <- countsVbFinalNonzeroQCBinsWIS[rownames(metaDataFinalNonzeroQC_HiSeq),]

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaDataFinalNonzeroQC_HiSeq_HMS <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metaDataFinalNonzeroQC_HiSeq_BCM <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metaDataFinalNonzeroQC_HiSeq_MDA <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metaDataFinalNonzeroQC_HiSeq_WashU <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metaDataFinalNonzeroQC_HiSeq_Broad_WGS <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metaDataFinalNonzeroQC_HiSeq_UNC <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metaDataFinalNonzeroQC_HiSeq_CMS <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

# WGS
metaDataFinalNonzeroQC_HiSeq_WGS <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA
metaDataFinalNonzeroQC_HiSeq_RNA <- metaDataFinalNonzeroQC_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset Full count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
countsVbFinalNonzeroQC_HiSeq_HMS <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_HMS),]
countsVbFinalNonzeroQC_HiSeq_BCM <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_BCM),]
countsVbFinalNonzeroQC_HiSeq_MDA <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_MDA),]
countsVbFinalNonzeroQC_HiSeq_WashU <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WashU),]
countsVbFinalNonzeroQC_HiSeq_Broad_WGS <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_Broad_WGS),]
countsVbFinalNonzeroQC_HiSeq_WGS <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
countsVbFinalNonzeroQC_HiSeq_UNC <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_UNC),]
countsVbFinalNonzeroQC_HiSeq_CMS <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_CMS),]
countsVbFinalNonzeroQC_HiSeq_RNA <- countsVbFinalNonzeroQC_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_RNA),]

#--------------------Subset WIS count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
countsVbFinalNonzeroQCWIS_HiSeq_HMS <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_HMS),]
countsVbFinalNonzeroQCWIS_HiSeq_BCM <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_BCM),]
countsVbFinalNonzeroQCWIS_HiSeq_MDA <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_MDA),]
countsVbFinalNonzeroQCWIS_HiSeq_WashU <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WashU),]
countsVbFinalNonzeroQCWIS_HiSeq_Broad_WGS <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_Broad_WGS),]
countsVbFinalNonzeroQCWIS_HiSeq_WGS <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
countsVbFinalNonzeroQCWIS_HiSeq_UNC <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_UNC),]
countsVbFinalNonzeroQCWIS_HiSeq_CMS <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_CMS),]
countsVbFinalNonzeroQCWIS_HiSeq_RNA <- countsVbFinalNonzeroQCWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_RNA),]

#--------------------Subset Bins count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
countsVbFinalNonzeroQCBins_HiSeq_HMS <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_HMS),]
countsVbFinalNonzeroQCBins_HiSeq_BCM <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_BCM),]
countsVbFinalNonzeroQCBins_HiSeq_MDA <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_MDA),]
countsVbFinalNonzeroQCBins_HiSeq_WashU <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WashU),]
countsVbFinalNonzeroQCBins_HiSeq_Broad_WGS <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_Broad_WGS),]
countsVbFinalNonzeroQCBins_HiSeq_WGS <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
countsVbFinalNonzeroQCBins_HiSeq_UNC <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_UNC),]
countsVbFinalNonzeroQCBins_HiSeq_CMS <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_CMS),]
countsVbFinalNonzeroQCBins_HiSeq_RNA <- countsVbFinalNonzeroQCBins_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_RNA),]

#--------------------Subset BinsWIS count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
countsVbFinalNonzeroQCBinsWIS_HiSeq_HMS <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_HMS),]
countsVbFinalNonzeroQCBinsWIS_HiSeq_BCM <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_BCM),]
countsVbFinalNonzeroQCBinsWIS_HiSeq_MDA <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_MDA),]
countsVbFinalNonzeroQCBinsWIS_HiSeq_WashU <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WashU),]
countsVbFinalNonzeroQCBinsWIS_HiSeq_Broad_WGS <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_Broad_WGS),]
countsVbFinalNonzeroQCBinsWIS_HiSeq_WGS <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
countsVbFinalNonzeroQCBinsWIS_HiSeq_UNC <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_UNC),]
countsVbFinalNonzeroQCBinsWIS_HiSeq_CMS <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_CMS),]
countsVbFinalNonzeroQCBinsWIS_HiSeq_RNA <- countsVbFinalNonzeroQCBinsWIS_HiSeq[rownames(metaDataFinalNonzeroQC_HiSeq_RNA),]

#--------------------Save data for ML--------------------#
save(# Subset raw count data
  countsVbFinalNonzeroQC_HiSeq_HMS,
  countsVbFinalNonzeroQC_HiSeq_BCM,
  countsVbFinalNonzeroQC_HiSeq_MDA,
  countsVbFinalNonzeroQC_HiSeq_WashU,
  countsVbFinalNonzeroQC_HiSeq_Broad_WGS,
  countsVbFinalNonzeroQC_HiSeq_UNC,
  countsVbFinalNonzeroQC_HiSeq_CMS,
  
  countsVbFinalNonzeroQC_HiSeq_WGS,
  countsVbFinalNonzeroQC_HiSeq_RNA,
  
  # Subset raw WIS count data
  countsVbFinalNonzeroQCWIS_HiSeq_HMS,
  countsVbFinalNonzeroQCWIS_HiSeq_BCM,
  countsVbFinalNonzeroQCWIS_HiSeq_MDA,
  countsVbFinalNonzeroQCWIS_HiSeq_WashU,
  countsVbFinalNonzeroQCWIS_HiSeq_Broad_WGS,
  countsVbFinalNonzeroQCWIS_HiSeq_UNC,
  countsVbFinalNonzeroQCWIS_HiSeq_CMS,
  
  countsVbFinalNonzeroQCWIS_HiSeq_WGS,
  countsVbFinalNonzeroQCWIS_HiSeq_RNA,
  
  # Subset raw Bins count data
  countsVbFinalNonzeroQCBins_HiSeq_HMS,
  countsVbFinalNonzeroQCBins_HiSeq_BCM,
  countsVbFinalNonzeroQCBins_HiSeq_MDA,
  countsVbFinalNonzeroQCBins_HiSeq_WashU,
  countsVbFinalNonzeroQCBins_HiSeq_Broad_WGS,
  countsVbFinalNonzeroQCBins_HiSeq_UNC,
  countsVbFinalNonzeroQCBins_HiSeq_CMS,
  
  countsVbFinalNonzeroQCBins_HiSeq_WGS,
  countsVbFinalNonzeroQCBins_HiSeq_RNA,
  
  # Subset raw BinsWIS count data
  countsVbFinalNonzeroQCBinsWIS_HiSeq_HMS,
  countsVbFinalNonzeroQCBinsWIS_HiSeq_BCM,
  countsVbFinalNonzeroQCBinsWIS_HiSeq_MDA,
  countsVbFinalNonzeroQCBinsWIS_HiSeq_WashU,
  countsVbFinalNonzeroQCBinsWIS_HiSeq_Broad_WGS,
  countsVbFinalNonzeroQCBinsWIS_HiSeq_UNC,
  countsVbFinalNonzeroQCBinsWIS_HiSeq_CMS,
  
  countsVbFinalNonzeroQCBinsWIS_HiSeq_WGS,
  countsVbFinalNonzeroQCBinsWIS_HiSeq_RNA,
  
  # Subset metadata
  metaDataFinalNonzeroQC_HiSeq_HMS,
  metaDataFinalNonzeroQC_HiSeq_BCM,
  metaDataFinalNonzeroQC_HiSeq_MDA,
  metaDataFinalNonzeroQC_HiSeq_WashU,
  metaDataFinalNonzeroQC_HiSeq_Broad_WGS,
  metaDataFinalNonzeroQC_HiSeq_UNC,
  metaDataFinalNonzeroQC_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metaDataFinalNonzeroQC_HiSeq_WGS,
  metaDataFinalNonzeroQC_HiSeq_RNA,
  
  # Full metadata
  metaDataFinalNonzeroQC_HiSeq,
  file = paste0("Interim_data/",plotPrefix,"/data_per_center_tcga_full_wis_bins_features_subset_20Feb24.RData"))

#----------------------------------------------------#
# Incorporate stage information
#----------------------------------------------------#
library(forcats)

# Remove "Not available" or minority path stages
metaDataFinalNonzeroQCPath <- droplevels(metaDataFinalNonzeroQC[! (metaDataFinalNonzeroQC$pathologic_stage_label == "Not available" | # n = 5950
                                                                     metaDataFinalNonzeroQC$pathologic_stage_label == "I or II NOS" | # n = 24
                                                                     metaDataFinalNonzeroQC$pathologic_stage_label == "Stage 0" | # n = 7
                                                                     metaDataFinalNonzeroQC$pathologic_stage_label == "Stage IS" | # n = 50
                                                                     metaDataFinalNonzeroQC$pathologic_stage_label == "Stage Tis" | # n = 1
                                                                     metaDataFinalNonzeroQC$pathologic_stage_label == "Stage X"),]) # n = 18

tumorStageVector <- metaDataFinalNonzeroQCPath$pathologic_stage_label
tumorStageVectorBinned <- fct_collapse(tumorStageVector,
                                       Stage1 = c("Stage I", "Stage IA", "Stage IB"),
                                       Stage2 = c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"),
                                       Stage3 = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),
                                       Stage4 = c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"))

# Add back to metadata
metaDataFinalNonzeroQCPath$pathologic_stage_label_binned <- tumorStageVectorBinned

## Identify early stage BDN samples
metaDataFinalNonzeroQCPath_BDN_EarlyStage <- metaDataFinalNonzeroQCPath %>%
  filter(sample_type == "Blood Derived Normal", 
         pathologic_stage_label_binned %in% c("Stage1","Stage2")) %>%
  droplevels()

#----------------------------------------------------#
# Identify blood samples without mutations
#----------------------------------------------------#

load("Input_data/metadataQC_noMutUpdatedMay172019.RData",verbose = TRUE)
# metadataSamplesAllQCCGC_noMutGuardant
# metadataSamplesAllQCCGC_noMutFoundationOne

metaDataFinalNonzeroQC_noMutGuardant <- metaDataFinalNonzeroQC %>%
  filter(knightlabID %in% rownames(metadataSamplesAllQCCGC_noMutGuardant),
         sample_type == "Blood Derived Normal") %>%
  droplevels()
dim(metaDataFinalNonzeroQC_noMutGuardant) # 472  41 (original was n=481)

metaDataFinalNonzeroQC_noMutFoundation <- metaDataFinalNonzeroQC %>%
  filter(knightlabID %in% rownames(metadataSamplesAllQCCGC_noMutFoundationOne),
         sample_type == "Blood Derived Normal") %>%
  droplevels()
dim(metaDataFinalNonzeroQC_noMutFoundation) # 281  41 (original was n=284)

save(metaDataFinalNonzeroQCPath,
     metaDataFinalNonzeroQCPath_BDN_EarlyStage,
     metaDataFinalNonzeroQC_noMutGuardant,
     metaDataFinalNonzeroQC_noMutFoundation,
     file = paste0("Interim_data/",plotPrefix,"/metadata_stage_noMut_subsets_20Feb24.RData"))

#----------------------------------------------------#
# Plot PCAs
#----------------------------------------------------#

pcaPlotting <- function(pcaObject,pcChoices, dataLabels, factorString, titleString){
  require(ggbiplot)
  require(ggsci)
  theme_update(plot.title = element_text(hjust = 0.5))
  g <- ggbiplot(pcaObject,pcChoices, obs.scale = 1, var.scale = 1,
                groups = dataLabels, ellipse = TRUE,
                alpha = 0.2,
                circle = TRUE,var.axes=FALSE) + 
    scale_color_nejm(name = factorString) +
    theme_bw() + 
    #theme(legend.direction = "horizontal", legend.position = "top") +
    ggtitle(titleString) + theme(plot.title = element_text(hjust = 0.5))
  
  print(g)
}

#-------------------Full-------------------#
pcaVoom <- pcaPlotting(pcaObject = prcomp(voomDataGenusKrakenQCFilt),
                       pcChoices = c(1,2),
                       dataLabels = metaDataFinalNonzeroQC$data_submitting_center_label,
                       factorString = "Sequencing Center",
                       titleString = "PCA w/o Batch Correction")
ggsave(plot = pcaVoom, 
       filename = paste0("Figures/",plotPrefix,"/pcaVoom_Full_Center.jpeg"),
       width = 16.2, height = 5.29, units = "in", dpi = "retina")

pcaVSNM_ST <- pcaPlotting(pcaObject = prcomp(vsnmDataGenusKrakenQCFilt),
                          pcChoices = c(1,2),
                          dataLabels = metaDataFinalNonzeroQC$data_submitting_center_label,
                          factorString = "Sequencing Center",
                          titleString = "PCA w/ SNM Correction\n(Target: Sample Type)")
ggsave(plot = pcaVSNM_ST, 
       filename = paste0("Figures/",plotPrefix,"/pcaVSNM_Full_Center.jpeg"),
       width = 16.2, height = 5.29, units = "in", dpi = "retina")

#-------------------WIS-------------------#
pcaVoom_WIS <- pcaPlotting(pcaObject = prcomp(voomDataGenusKrakenQCFiltWIS),
                           pcChoices = c(1,2),
                           dataLabels = metaDataFinalNonzeroQC$data_submitting_center_label,
                           factorString = "Sequencing Center",
                           titleString = "PCA w/o Batch Correction")
ggsave(plot = pcaVoom_WIS, 
       filename = paste0("Figures/",plotPrefix,"/pcaVoom_WIS_Center.jpeg"),
       width = 16.2, height = 5.29, units = "in", dpi = "retina")

pcaVSNM_WIS_ST <- pcaPlotting(pcaObject = prcomp(vsnmDataGenusKrakenQCFiltWIS),
                              pcChoices = c(1,2),
                              dataLabels = metaDataFinalNonzeroQC$data_submitting_center_label,
                              factorString = "Sequencing Center",
                              titleString = "PCA w/ SNM Correction\n(Target: Sample Type)")
ggsave(plot = pcaVSNM_WIS_ST, 
       filename = paste0("Figures/",plotPrefix,"/pcaVSNM_WIS_Center.jpeg"),
       width = 16.2, height = 5.29, units = "in", dpi = "retina")

#-------------------Bins-------------------#
pcaVoom_Bins <- pcaPlotting(pcaObject = prcomp(voomDataGenusKrakenQCFiltBins),
                            pcChoices = c(1,2),
                            dataLabels = metaDataFinalNonzeroQC$data_submitting_center_label,
                            factorString = "Sequencing Center",
                            titleString = "PCA w/o Batch Correction")
ggsave(plot = pcaVoom_Bins, 
       filename = paste0("Figures/",plotPrefix,"/pcaVoom_Bins_Center.jpeg"),
       width = 16.2, height = 5.29, units = "in", dpi = "retina")

pcaVSNM_Bins_ST <- pcaPlotting(pcaObject = prcomp(vsnmDataGenusKrakenQCFiltBins),
                               pcChoices = c(1,2),
                               dataLabels = metaDataFinalNonzeroQC$data_submitting_center_label,
                               factorString = "Sequencing Center",
                               titleString = "PCA w/ SNM Correction\n(Target: Sample Type)")
ggsave(plot = pcaVSNM_Bins_ST, 
       filename = paste0("Figures/",plotPrefix,"/pcaVSNM_Bins_Center.jpeg"),
       width = 16.2, height = 5.29, units = "in", dpi = "retina")

#----------------------------------------------------#
# Calculate PVCAs
#----------------------------------------------------#
require(lme4)

source("00-functions.R") # for PVCA() function
pct_threshold <- 0.8
metaPVCAExtendedFiltered <- metaDataFinalNonzeroQC[,c("sample_type",
                                                      "disease_type",
                                                      "data_submitting_center_label",
                                                      "cgc_platform",
                                                      "experimental_strategy", 
                                                      "tissue_source_site_label",
                                                      "portion_is_ffpe")]

#-------------------Full-------------------#
pvcaRaw_Full <- PVCA(counts = t(countsVbFinalNonzeroQC), 
                     meta = metaPVCAExtendedFiltered, 
                     threshold = pct_threshold,
                     inter = FALSE)
pvcaRaw_Full

pvcaVoom_Full <- PVCA(counts = t(voomDataGenusKrakenQCFilt),
                      meta = metaPVCAExtendedFiltered,
                      threshold = pct_threshold,
                      inter = FALSE)
pvcaVoom_Full

pvcaVSNM_Full <- PVCA(counts = t(vsnmDataGenusKrakenQCFilt), 
                      meta = metaPVCAExtendedFiltered,
                      threshold = pct_threshold,
                      inter = FALSE)
pvcaVSNM_Full

#-------------------WIS-------------------#
pvcaRaw_WIS <- PVCA(counts = t(countsVbFinalNonzeroQCWIS), 
                    meta = metaPVCAExtendedFiltered, 
                    threshold = pct_threshold,
                    inter = FALSE)
pvcaRaw_WIS

pvcaVoom_WIS <- PVCA(counts = t(voomDataGenusKrakenQCFiltWIS),
                     meta = metaPVCAExtendedFiltered,
                     threshold = pct_threshold,
                     inter = FALSE)
pvcaVoom_WIS

pvcaVSNM_WIS <- PVCA(counts = t(vsnmDataGenusKrakenQCFiltWIS), 
                     meta = metaPVCAExtendedFiltered,
                     threshold = pct_threshold,
                     inter = FALSE)
pvcaVSNM_WIS
save(pvcaVSNM_WIS, file = paste0("Interim_data/",plotPrefix,"/pvcaVSNM_WIS.RData"))

#----------------------------------------------------#
# Correlation of counts
#----------------------------------------------------#

sampleIDOverlap <- intersect(metaDataFinalNonzeroQC$knightlabID,
                             rownames(vbDataBarnDFReconciledQC))
featOverlap <- intersect(colnames(countsVbFinalNonzeroQC),
                         colnames(vbDataBarnDFReconciledQC))
countsVbFinalNonzeroQCFormTmp <- countsVbFinalNonzeroQC
rownames(countsVbFinalNonzeroQCFormTmp) <- metaDataFinalNonzeroQC$knightlabID
countsVbFinalNonzeroQCForm <- countsVbFinalNonzeroQCFormTmp[sampleIDOverlap,featOverlap]
# Realign orig data
vbDataBarnDFReconciledQC_Matched <- vbDataBarnDFReconciledQC[sampleIDOverlap,featOverlap]

dim(countsVbFinalNonzeroQCForm) # 15180  1885
dim(vbDataBarnDFReconciledQC_Matched) # 15180  1885

#--------------------Corr--------------------#
require(ggpubr)
## By sample counts
cor.test(log10(rowSums(vbDataBarnDFReconciledQC_Matched)),
         log10(rowSums(countsVbFinalNonzeroQCForm)),
         method = "spearman") # 0.6613132

# Sanity check
all(rownames(countsVbFinalNonzeroQCForm) == metaDataFinalNonzeroQC$knightlabID) # TRUE
dfSampleCounts <- data.frame(Original = log10(rowSums(vbDataBarnDFReconciledQC_Matched)),
                             HPRC = log10(rowSums(countsVbFinalNonzeroQCForm)),
                             experimental_strategy = metaDataFinalNonzeroQC$experimental_strategy)

dfSampleCounts %>%
  ggscatter(x = "Original",
            y = "HPRC",
            xlab = "Original TCGA sample counts (log10)",
            ylab = "HPRC-depleted TCGA sample counts (log10)",
            # add = "reg.line",
            # add.params = list(color = "blue"),
            # conf.int = TRUE,
            alpha = 0.3,
            # facet.by = "experimental_strategy",
            cor.coef = FALSE,
            cor.method = "spearman") +
  theme(aspect.ratio=1) +
  stat_cor(method = "spearman",cor.coef.name = "rho")
ggsave(filename = paste0("Figures/",plotPrefix,"/corr_by_sample.jpeg"),
       width = 5, height = 5, units = "in", dpi = "retina")

## By microbe
# Remove zero-sum microbes
cor.test(log10(colSums(vbDataBarnDFReconciledQC_Matched)),
         log10(colSums(countsVbFinalNonzeroQCForm)),
         method = "spearman") # 0.7741788

dfMicrobeCounts <- data.frame(Original = log10(colSums(vbDataBarnDFReconciledQC_Matched)),
                              HPRC = log10(colSums(countsVbFinalNonzeroQCForm))) %>%
  rownames_to_column("Taxa") %>%
  mutate(Taxa = gsub("^k__.+\\.g__","",Taxa)) %>%
  mutate(WIS = ifelse(Taxa %in% wisGeneraBactKnownUniqueRev,TRUE,FALSE)) %>%
  column_to_rownames("Taxa")

dfMicrobeCounts %>%
  ggscatter(x = "Original",
            y = "HPRC",
            color = "WIS",
            palette = c("darkgrey", "red"),
            xlab = "Original TCGA microbial counts (log10)",
            ylab = "HPRC-depleted TCGA microbial counts (log10)",
            # add = "reg.line",
            # add.params = list(color = "blue"),
            # conf.int = TRUE,
            alpha = 0.2,
            # facet.by = "experimental_strategy",
            cor.coef = FALSE) +
  theme(aspect.ratio=1) +
  stat_cor(method = "spearman",cor.coef.name = "rho")
ggsave(filename = paste0("Figures/",plotPrefix,"/corr_by_microbe.jpeg"),
       width = 5, height = 5, units = "in", dpi = "retina")

dfMicrobeCounts %>%
  filter(WIS == TRUE) %>%
  ggscatter(x = "Original",
            y = "HPRC",
            xlab = "Original WIS TCGA microbial counts (log10)",
            ylab = "HPRC-depleted WIS TCGA microbial counts (log10)",
            # add = "reg.line",
            # add.params = list(color = "blue"),
            # conf.int = TRUE,
            alpha = 0.2,
            # facet.by = "experimental_strategy",
            cor.coef = FALSE) +
  theme(aspect.ratio=1) +
  stat_cor(method = "spearman",cor.coef.name = "rho")
ggsave(filename = paste0("Figures/",plotPrefix,"/corr_by_microbe_WISonly.jpeg"),
       width = 5, height = 5, units = "in", dpi = "retina")

#----------------------------------------------------#
# Plot ML results
#----------------------------------------------------#
require(tidyr)
require(pheatmap)

abbreviationsTCGA_Allcancer <- read.csv("Input_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

mlPerfAll <- read.csv("Supporting_scripts/01v2-ML-Pangenome/perfML_Pangenome_Full_WIS_ALL_13Feb24.csv", 
                      stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=aucroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

# Format variable names
mlPerfAll$datasetName[mlPerfAll$datasetName == "vsnmDataGenusKrakenQCFilt"] <- "Full"
mlPerfAll$datasetName[mlPerfAll$datasetName == "vsnmDataGenusKrakenQCFiltWIS"] <- "WIS"
mlPerfAll$datasetName <- factor(mlPerfAll$datasetName, levels = c("Full", "WIS"))
table(mlPerfAll$datasetName)
table(mlPerfAll$metadataName)

## Heatmaps
heatmapPlot <- function(data=mlPerfAll, st, dt,meta=NULL,fontSize=16){
  if(is.null(meta)){
    tmp <- data %>%
      filter(sampleType == st,
             datasetName == dt) %>%
      select(abbrev, AUROC, AUPR) %>%
      distinct(abbrev, .keep_all=TRUE) %>%
      arrange(abbrev) %>%
      column_to_rownames("abbrev")
  } else{
    tmp <- data %>%
      filter(sampleType == st,
             datasetName == dt,
             metadataName %in% meta) %>%
      select(abbrev, AUROC, AUPR) %>%
      distinct(abbrev, .keep_all=TRUE) %>%
      arrange(abbrev) %>%
      column_to_rownames("abbrev")
  }
  tmp %>% t() %>%
    pheatmap(color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(50),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             display_numbers = TRUE,
             number_format = "%.2f",
             na_col = "white",
             filename = paste0("Figures/",plotPrefix,"/heatmap_",
                               dt,"_",st,"__",paste(meta,sep="_"),".jpeg"),
             fontsize = fontSize,
             number_color = "black",
             width = 20,
             height = 2)
  tmp %>% summarise(meanAUROC = mean(AUROC),
                    meanAUPR = mean(AUPR)) %>% print()
}

#----------------Full----------------#
## PT
heatmapPlot("Primary Tumor", "Full", "metaDataFinalNonzeroQC")
## BDN
heatmapPlot("Blood Derived Normal", "Full", "metaDataFinalNonzeroQC",18)
## STN
heatmapPlot("Primary Tumor vs Solid Tissue Normal", "Full", "metaDataFinalNonzeroQC",18)

#----------------WIS----------------#
## PT
heatmapPlot("Primary Tumor", "WIS", "metaDataFinalNonzeroQC")
## BDN
heatmapPlot("Blood Derived Normal", "WIS", "metaDataFinalNonzeroQC",18)
## STN
heatmapPlot("Primary Tumor vs Solid Tissue Normal", "WIS", "metaDataFinalNonzeroQC",18)

#----------------Full: Early stage and noMut BDN----------------#
## Early stage
heatmapPlot("Blood Derived Normal", "Full", "metaDataFinalNonzeroQCPath_BDN_EarlyStage",18)
## noMut Guardant
heatmapPlot("Blood Derived Normal", "Full", "metaDataFinalNonzeroQC_noMutGuardant",18)
## noMut Foundation
heatmapPlot("Blood Derived Normal", "Full", "metaDataFinalNonzeroQC_noMutFoundation",18)

#----------------WIS: Early stage and noMut BDN----------------#
## Early stage
heatmapPlot("Blood Derived Normal", "WIS", "metaDataFinalNonzeroQCPath_BDN_EarlyStage",18)
## noMut Guardant
heatmapPlot("Blood Derived Normal", "WIS", "metaDataFinalNonzeroQC_noMutGuardant",18)
## noMut Foundation
heatmapPlot("Blood Derived Normal", "WIS", "metaDataFinalNonzeroQC_noMutFoundation",18)

#----------------------------------------------------#
# Plot ML results: WGS vs RNA
#----------------------------------------------------#
require(tidyr)
require(pheatmap)

abbreviationsTCGA_Allcancer <- read.csv("Input_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

mlPerfAll_WGSvsRNA <- read.csv("Supporting_scripts/01v3-ML-Pangenome/perfML_Pangenome_Full_WIS_WGSvsRNA_ALL_13Feb24.csv", 
                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=aucroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

# Format variable names
mlPerfAll_WGSvsRNA$datasetName[mlPerfAll_WGSvsRNA$datasetName == "vsnmDataGenusKrakenQCFilt"] <- "Full"
mlPerfAll_WGSvsRNA$datasetName[mlPerfAll_WGSvsRNA$datasetName == "vsnmDataGenusKrakenQCFiltWIS"] <- "WIS"
mlPerfAll_WGSvsRNA$datasetName <- factor(mlPerfAll_WGSvsRNA$datasetName, levels = c("Full", "WIS"))
table(mlPerfAll_WGSvsRNA$datasetName)
table(mlPerfAll_WGSvsRNA$metadataName)

#----------------Full----------------#
## PT, WGS
heatmapPlot(mlPerfAll_WGSvsRNA, "Primary Tumor", "Full", "metaDataFinalNonzeroQC_WGS")
## PT, RNA
heatmapPlot(mlPerfAll_WGSvsRNA, "Primary Tumor", "Full", "metaDataFinalNonzeroQC_RNA")

## STN, WGS
heatmapPlot(mlPerfAll_WGSvsRNA, "Primary Tumor vs Solid Tissue Normal", "Full", "metaDataFinalNonzeroQC_WGS")
## STN, RNA
heatmapPlot(mlPerfAll_WGSvsRNA, "Primary Tumor vs Solid Tissue Normal", "Full", "metaDataFinalNonzeroQC_RNA")


## BDN
heatmapPlot(mlPerfAll_WGSvsRNA, "Blood Derived Normal", "Full", "metaDataFinalNonzeroQC_WGS",18)
## STN
heatmapPlot(mlPerfAll_WGSvsRNA, "Primary Tumor vs Solid Tissue Normal", "Full", "metaDataFinalNonzeroQC",18)

#----------------WIS----------------#
## PT
heatmapPlot("Primary Tumor", "WIS", "metaDataFinalNonzeroQC")
## BDN
heatmapPlot("Blood Derived Normal", "WIS", "metaDataFinalNonzeroQC",18)
## STN
heatmapPlot("Primary Tumor vs Solid Tissue Normal", "WIS", "metaDataFinalNonzeroQC",18)



#----------------------------------------------------#
# Original ecological validation analyses
#----------------------------------------------------#
require(ggpubr)
require(ggsci)

load(paste0("Interim_data/",plotPrefix,"/data_vsnm_tcga_full_and_wis_features_subset_12Feb24.RData"),
     verbose = TRUE)
# countsVbFinalNonzeroQC
# voomDataGenusKrakenQCFilt
# vsnmDataGenusKrakenQCFilt
# countsVbFinalNonzeroQCWIS
# voomDataGenusKrakenQCFiltWIS
# vsnmDataGenusKrakenQCFiltWIS
# metaDataFinalNonzeroQC

load("Input_data/cgcAPIMetadataJoined.RData", verbose = TRUE)
# metadataSamplesAllQCCGC
load("Input_data/cgcMetadataKrakenProj.RData", verbose = TRUE)
# cgcMetadataKrakenProj

# # See: https://cran.r-project.org/src/contrib/Archive/bigrquery/
# install.packages("https://cran.r-project.org/src/contrib/Archive/bigrquery/bigrquery_1.2.0.tar.gz",
#                  repos = NULL, type = "source")
require(bigrquery)

#------------- HPV status across CESC -------------#

clinical_table = "[isb-cgc:tcga_201607_beta.Clinical_data]"
cloud_project_workshop = "hybrid-coyote-219120"
sqlQuery = paste("SELECT ParticipantBarcode, Study, hpv_calls, hpv_status ", 
                 "FROM ", clinical_table, sep="")
sqlQuery
hpv_table = query_exec(sqlQuery,project = cloud_project_workshop)
# save(hpv_table,
#      file = "Interim_data/Pangenome/hpv_table_12Feb24.RData")

isbcgcHPV <- hpv_table
hpvPancancerMeta <- left_join(metaDataFinalNonzeroQC, 
                              isbcgcHPV, 
                              by = c("tcga_case_id" = "ParticipantBarcode") )
hpvPancancerData <- data.frame(HPV = countsVbFinalNonzeroQC[rownames(metaDataFinalNonzeroQC),
                                                            "k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
hpvPancancerCombined <- droplevels(cbind(hpvPancancerMeta, hpvPancancerData))

interactVec <- as.character(interaction(hpvPancancerCombined$investigation, 
                                        hpvPancancerCombined$hpv_status, sep = " "))
interactVec[which(is.na(interactVec))] <- as.character(hpvPancancerCombined$investigation[which(is.na(interactVec))])
interactVec <- factor(interactVec)
hpvPancancerCombined$hpvInteract <- interactVec

hpvCervicalCancerComparisons <- list( c("TCGA-CESC Positive", "TCGA-CESC Negative"))
hpvPancancerCombined %>%
  filter((sample_type %in% c("Blood Derived Normal", "Primary Tumor")) &
           !(hpv_status %in% c("Indeterminate")) &
           (disease_type == "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma")) %>%
  ggboxplot(x = "hpvInteract", y = "HPV", 
            # color = "sample_type",
            add = "jitter",
            facet.by = "sample_type",
            # palette = pal_nejm(),
            xlab = "Clinical HPV Status", ylab = "SNM Normalized Abundance",
            # ylim = c(-11, 20),
            title = "Pancancer Comparison of Alphapapillomavirus Genus Abundance in Cervical Cancer",
            # legend = "right",
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels=c("TCGA-CESC Negative" = "Negative",
                            "TCGA-CESC Positive" = "Positive")) +
  scale_color_nejm() +
  rotate_x_text(angle = 30) +
  scale_y_log10() +
  stat_compare_means(comparisons = hpvCervicalCancerComparisons, 
                     label = "p.signif", 
                     method.args = list(alternative = "greater")) #-> p # Add pairwise comparisons p-value

#------------- LCV status in stomach cancers -------------#

stadMasterPatientTable <- read.csv(file = "Input_data/STAD_Master_Patient_Table_20140207.csv", 
                                   stringsAsFactors = FALSE)

stadPancancerMeta <- left_join(metaDataFinalNonzeroQC,
                               stadMasterPatientTable, 
                               by = c("tcga_case_id" = "TCGA.barcode") )
rownames(stadPancancerMeta) <- rownames(metaDataFinalNonzeroQC)
lcvPancancerData <- data.frame(LCV = vsnmDataGenusKrakenQCFilt[rownames(stadPancancerMeta),"k__Viruses.o__Herpesvirales.f__Herpesviridae.g__Lymphocryptovirus"],
                               HPylori = vsnmDataGenusKrakenQCFilt[rownames(stadPancancerMeta),"k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae.g__Helicobacter"])
lcvPancancerCombined <- droplevels(cbind(stadPancancerMeta, lcvPancancerData))

lcvComparisons <- list( c("EBV", "CIN"), c("EBV", "GS"), c("EBV","MSI") )
lcvPancancerCombined %>%
  filter((sample_type %in% c("Blood Derived Normal", "Solid Tissue Normal", "Primary Tumor")) &
           !(is.na(Molecular.Subtype)) &
           (disease_type == "Stomach Adenocarcinoma")) %>%
  ggboxplot(x = "Molecular.Subtype", y = "LCV", 
            # color = "sample_type",
            add = "jitter",
            facet.by = "sample_type",
            # palette = "lancet",
            ylim = c(-5, 22),
            xlab = "STAD Molecular subtype (The Cancer Genome Atlas Research Network, 2014. Nature)", ylab = "SNM Normalized Abundance", 
            title = "Pancancer Comparison of Lymphocryptovirus Genus Abundance in Stomach Adenocarcinoma",
            # legend = "right",
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  scale_color_nejm() +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = lcvComparisons, label = "p.signif", method = "wilcox.test") # -> p # Add pairwise comparisons p-value

ggsave(p, filename = "EBV in STAD.png", path = "./Clinical Validation Plots",dpi = "retina", units = "in",
       height = 4, width = 7)

#------------- HPV status in HNSC -------------#
# Install from unzipped zip file of Github repo
require(TCGAmutations)
hnscClinicalMetadata <- as.data.frame(getClinicalData(tcga_load(study = "HNSC"))) #

# Data alignment
hnscClinicalCols2Keep <- c(
  "Tumor_Sample_Barcode",
  "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
  "hpv_status_by_p16_testing",
  "hpv_status_by_ish_testing"
)

hnscMetadataQCCGC <- metaDataFinalNonzeroQC[metaDataFinalNonzeroQC$disease_type == "Head and Neck Squamous Cell Carcinoma",]
hnscMetadataQCCGC$case_uuid <- toupper(hnscMetadataQCCGC$cgc_case_uuid)
hnscClinicalMetadata$bcr_patient_uuid <- toupper(hnscClinicalMetadata$bcr_patient_uuid)

hnscMetadataQCCGCClinical <- left_join(hnscMetadataQCCGC, 
                                       hnscClinicalMetadata[,hnscClinicalCols2Keep], 
                                       by = c("cgc_case_uuid" = "bcr_patient_uuid"))
rownames(hnscMetadataQCCGCClinical) <- rownames(hnscMetadataQCCGC)

# Subset data
hnscMetadataQCCGCClinical_HPVp16 <- droplevels(hnscMetadataQCCGCClinical[!is.na(hnscMetadataQCCGCClinical$hpv_status_by_p16_testing),])
hnscMetadataQCCGCClinical_HPVish <- droplevels(hnscMetadataQCCGCClinical[!is.na(hnscMetadataQCCGCClinical$hpv_status_by_ish_testing),])

hnscHPVp16Data <- data.frame(HPV = vsnmDataGenusKrakenQCFilt[rownames(hnscMetadataQCCGCClinical_HPVp16),"k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
hnscHPVishData <- data.frame(HPV = vsnmDataGenusKrakenQCFilt[rownames(hnscMetadataQCCGCClinical_HPVish),"k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
hnscHPVp16Combined <- cbind(hnscMetadataQCCGCClinical_HPVp16, hnscHPVp16Data)
hnscHPVishCombined <- cbind(hnscMetadataQCCGCClinical_HPVish, hnscHPVishData)

testType <- factor(c(rep("p16 Testing",dim(hnscHPVp16Combined)[1]), rep("ISH Testing", dim(hnscHPVishCombined)[1])))
testValue <- factor(c(as.character(hnscHPVp16Combined$hpv_status_by_p16_testing), as.character(hnscHPVishCombined$hpv_status_by_ish_testing)))
hnscHPVbothCombined <- cbind(rbind(hnscHPVp16Combined,hnscHPVishCombined),testType, testValue)
hnscHPVcomparisons <- list( c("Negative", "Positive"))
hnscHPVbothCombined %>%
  filter(sample_type == "Primary Tumor") %>%
  filter(!is.na(hpv_status_by_p16_testing)) %>%
  filter(!is.na(hpv_status_by_ish_testing)) %>%
  ggboxplot(x = "testValue", y = "HPV", 
            # color = "testValue",
            facet.by = "testType",
            add = "jitter",
            # palette = "lancet",
            xlab = "Clinical Testing for HPV", ylab = "SNM Normalized Abundance", 
            ylim = c(-3, 18),
            title = "Comparison of Alphapapillomavirus Genus Abundance in\nHead and Neck Squamous Cell Carcinoma Primary Tumors",
            legend = "right",
            legend.title = "Clinical HPV Status",
            font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_nejm() +
  rotate_x_text(angle = 30) +
  stat_compare_means(comparisons = hnscHPVcomparisons, 
                     label = "p.signif",
                     method.args = list(alternative = "less"),
                     method = "t.test") # -> p # Add pairwise comparisons p-value

hnscHPVcomparisons <- list( c("Negative", "Positive"))
hnscHPVp16Combined %>%
  filter(sample_type == "Primary Tumor") %>%
  ggboxplot(x = "hpv_status_by_p16_testing", y = "HPV", 
            color = "hpv_status_by_p16_testing",
            add = "jitter",
            palette = "lancet",
            xlab = "Clinical Testing for HPV", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Alphapapillomavirus Genus Abundance in\nHead and Neck Squamous Cell Carcinoma Primary Tumors",
            legend = "right",
            legend.title = "HPV Status by\nP16 Testing",
            font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = hnscHPVcomparisons, 
                     method.args = list(alternative = "less"),
                     method = "t.test") # Add pairwise comparisons p-value

hnscHPVcomparisons <- list( c("Negative", "Positive"))
hnscHPVishCombined %>%
  filter(sample_type == "Primary Tumor") %>%
  ggboxplot(x = "hpv_status_by_ish_testing", y = "HPV", 
            color = "hpv_status_by_ish_testing",
            add = "jitter",
            palette = "lancet",
            xlab = "Clinical Testing for HPV", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Alphapapillomavirus Genus Abundance in\nHead and Neck Squamous Cell Carcinoma Primary Tumors",
            legend = "right",
            legend.title = "HPV Status by\nISH Testing",
            font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = hnscHPVcomparisons, 
                     method.args = list(alternative = "less"),
                     method = "t.test") # Add pairwise comparisons p-value

#------------- HBV/HCV status in LIHC -------------#
lihcClinicalMetadata <- as.data.frame(getClinicalData(tcga_load(study = "LIHC"))) %>%
  filter(history_hepato_carcinoma_risk_factors %in% 
           c("Hepatitis_B","Alcohol_consumption","Hepatitis_C")) %>%
  filter(bcr_patient_uuid != "6dcfc418-ab28-4365-95ea-c1ae254f2341")

# Data alignment
lihcClinicalCols2Keep <- c(
  "Tumor_Sample_Barcode",
  "bcr_patient_uuid", ## NB: This aligns with the case_uuid in the QIIME mapping file
  "history_hepato_carcinoma_risk_factors"
)

lihcMetadataQCCGC <- metaDataFinalNonzeroQC[metaDataFinalNonzeroQC$disease_type == "Liver Hepatocellular Carcinoma",] 
lihcMetadataQCCGC <- metaDataFinalNonzeroQC %>% 
  filter(disease_type == "Liver Hepatocellular Carcinoma") %>%
  filter(cgc_case_uuid != "D6486001-240A-455A-980C-E06C25C61FA5")
lihcMetadataQCCGC$case_uuid <- toupper(lihcMetadataQCCGC$cgc_case_uuid)
lihcClinicalMetadata$bcr_patient_uuid <- toupper(lihcClinicalMetadata$bcr_patient_uuid)

lihcMetadataQCCGCClinical <- left_join(lihcMetadataQCCGC, 
                                       lihcClinicalMetadata[,lihcClinicalCols2Keep], 
                                       by = c("case_uuid" = "bcr_patient_uuid"))
rownames(lihcMetadataQCCGCClinical) <- rownames(lihcMetadataQCCGC)

# Subset data
lihcMetadataQCCGCClinical_Riskfactors <- droplevels(lihcMetadataQCCGCClinical[!is.na(lihcMetadataQCCGCClinical$history_hepato_carcinoma_risk_factors),])

hbvHcvGenera <- c("k__Viruses.f__Flaviviridae.g__Hepacivirus", "k__Viruses.f__Hepadnaviridae.g__Orthohepadnavirus")
lihcHepData <- vsnmDataGenusKrakenQCFilt[rownames(lihcMetadataQCCGCClinical_Riskfactors),hbvHcvGenera]
lihcHepDataCombined <- droplevels(cbind(lihcMetadataQCCGCClinical_Riskfactors, lihcHepData))

# my_comparisons <- list( c("Primary Tumor", "Solid Tissue Normal"), c("Primary Tumor", "Blood Derived Normal"), c("Solid Tissue Normal", "Blood Derived Normal"))
lihcHepComparisons <- list( c("Hepatitis_B", "Hepatitis_C"), c("Hepatitis_B","Alcohol_consumption"), c("Hepatitis_C","Alcohol_consumption"))
lihcHepDataCombined %>%
  filter((history_hepato_carcinoma_risk_factors %in% c("Hepatitis_C", "Hepatitis_B","Alcohol_consumption")) &
           (sample_type %in% c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"))) %>%
  ggboxplot(x = "history_hepato_carcinoma_risk_factors", 
            y = "k__Viruses.f__Hepadnaviridae.g__Orthohepadnavirus", 
            # color = "sample_type",
            facet.by = "sample_type",
            palette = "lancet",
            add = c("jitter"), # add = "jitter", 
            # shape = "sample_type",
            xlab = "Clinically Assessed Patient History Risk Factors for Hepatocellular Carcinoma", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Orthohepadnavirus Genus Abundance in Liver Hepatocellular Carcinoma",
            legend = "right",
            legend.title = "Sample Type",
            font.label = list(size = 14, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("Alcohol_consumption" = "EtOH", 
                            "Hepatitis_B" = "Hep B",
                            "Hepatitis_C" = "Hep C")) +
  # rotate_x_text(angle = 30) +
  stat_compare_means(comparisons = lihcHepComparisons, label = "p.signif") # -> p # Add pairwise comparisons p-value

#-------------------------------#
# Fusobacterium GI-tract cancers vs others
#-------------------------------#
require(ggplot2)
require(ggpubr)
require(reshape2)
require(ggsci)
require(purrr)
require(ggrepel)

coadFusoMeta <- metaDataFinalNonzeroQC
coadFusoData <- data.frame(Fuso = vsnmDataGenusKrakenQCFilt[rownames(coadFusoMeta),"k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Fusobacteriaceae.g__Fusobacterium"])
coadFusoCombined <- cbind(coadFusoMeta, coadFusoData)

coadFusoCombined$GIBoolean <- factor(ifelse(coadFusoCombined$disease_type %in% c("Colon Adenocarcinoma",
                                                                                 "Rectum Adenocarcinoma",
                                                                                 "Cholangiocarcinoma",
                                                                                 "Liver Hepatocellular Carcinoma",
                                                                                 "Pancreatic Adenocarcinoma",
                                                                                 "Head and Neck Squamous Cell Carcinoma",
                                                                                 "Esophageal Carcinoma",
                                                                                 "Stomach Adenocarcinoma"),
                                            yes = "GI-Tract Cancer", no = "Non-GI-Tract Cancer"))

coadFusoGIComp <- list( c("GI-Tract Cancer", "Non-GI-Tract Cancer"))
coadFusoCombined %>%
  # filter(experimental_strategy %in% c("RNA-Seq")) %>%
  # filter(pathologic_stage_label %in% c("Stage IA","Stage IIA","Stage IIIA","Stage IVA")) %>%
  # filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  # filter(sample_type %in% c("Primary Tumor")) %>%
  # filter(disease_type %in% c("Colon Adenocarcinoma", "Stomach Adenocarcinoma", "Esophageal Carcinoma", "Rectum Adenocarcinoma")) %>%
  # filter(disease_type %in% c("Rectum Adenocarcinoma")) %>%
  ggboxplot(x = "GIBoolean", y = "Fuso", 
            # color = "GIBoolean",
            add = "jitter",
            add.params = list(alpha = 0.1),
            palette = "lancet",
            facet.by = "sample_type",
            nrow = 1,
            xlab = "Cancer Type", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Fusobacterium Genus Abundance",
            legend = "right",
            ylim = c(-1, 21),
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  rotate_x_text(angle = 30) +
  theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
  stat_summary(geom = "text", 
               fun.data = function(x){c(y = -0.5, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0.75)) +
  stat_compare_means(comparisons = coadFusoGIComp, 
                     label.y = c(20,20,20),
                     #size = 8,
                     label = "p.format") -> fusoPlotGI # Add pairwise comparisons p-value

coadFusoCombined %>%
  # filter(experimental_strategy %in% c("RNA-Seq")) %>%
  # filter(pathologic_stage_label %in% c("Stage IA","Stage IIA","Stage IIIA","Stage IVA")) %>%
  # filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  # filter(sample_type %in% c("Primary Tumor")) %>%
  # filter(disease_type %in% c("Colon Adenocarcinoma", "Stomach Adenocarcinoma", "Esophageal Carcinoma", "Rectum Adenocarcinoma")) %>%
  # filter(disease_type %in% c("Rectum Adenocarcinoma")) %>%
  ggboxplot(x = "GIBoolean", y = "Fuso", 
            # color = "GIBoolean",
            add = "jitter",
            add.params = list(alpha = 0.05),
            palette = "lancet",
            facet.by = "sample_type",
            nrow = 1,
            xlab = "Cancer Type", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Fusobacterium Genus Abundance",
            legend = "right",
            ylim = c(-1, 30),
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  rotate_x_text(angle = 30) +
  theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 14)) -> fusoPlotGI # Add pairwise comparisons p-value







#--------------------------------------------------------------------------------------------------------------------#
# Separate raw data into seq center-experimental strategy groups (to preclude needing batch correction)
#--------------------------------------------------------------------------------------------------------------------#
metadataSamplesAllQC %>% count(data_submitting_center_label, experimental_strategy)
metadataSamplesAllQC %>% count(platform) # HiSeq accounts for 91.27% of samples
metadataSamplesAllQC %>% count(platform, data_submitting_center_label)

# Subset metadata to Illumina HiSeq samples
metadataSamplesAll_HiSeq <- metadataSamplesAll %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()
metadataSamplesAllQC_HiSeq <- metadataSamplesAllQC %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()

# Subset count data to Illumina HiSeq samples
tcgaGenusKrakenAllFiltWIS_HiSeq <- tcgaGenusKrakenAllFiltWIS[rownames(metadataSamplesAll_HiSeq),]
tcgaGenusKrakenQCFiltWIS_HiSeq <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq),]

dim(tcgaGenusKrakenAllFiltWIS_HiSeq) # 16243 184
dim(tcgaGenusKrakenQCFiltWIS_HiSeq) # 16087 184

# save(tcgaGenusKrakenAllFiltWIS_HiSeq,
#      tcgaGenusKrakenQCFiltWIS_HiSeq,
#      metadataSamplesAll_HiSeq,
#      metadataSamplesAllQC_HiSeq,
#      file = "Interim_data/tcga-hiseq-wis-overlapping-data-and-metadata-subset-25July22.RData")

#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_HMS <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metadataSamplesAllQC_HiSeq_BCM <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Baylor College of Medicine") %>% droplevels()
metadataSamplesAllQC_HiSeq_MDA <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% droplevels()
metadataSamplesAllQC_HiSeq_WashU <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Washington University School of Medicine") %>% droplevels()
metadataSamplesAllQC_HiSeq_Broad_WGS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_UNC <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()
metadataSamplesAllQC_HiSeq_CMS <- metadataSamplesAllQC_HiSeq %>% filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% droplevels()
metadataSamplesAllQC_HiSeq_Broad_RNA <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_HMS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_HMS),]
tcgaGenusKrakenQCFiltWIS_BCM <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_BCM),]
tcgaGenusKrakenQCFiltWIS_MDA <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_MDA),]
tcgaGenusKrakenQCFiltWIS_WashU <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_WashU),]
tcgaGenusKrakenQCFiltWIS_Broad_WGS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
tcgaGenusKrakenQCFiltWIS_UNC <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_UNC),]
tcgaGenusKrakenQCFiltWIS_CMS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_CMS),]
tcgaGenusKrakenQCFiltWIS_Broad_RNA <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_Broad_RNA),]

#--------------------Subset metadata and count data by WGS (for multi-class classification)--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_WGS <- metadataSamplesAllQC_HiSeq %>% filter(experimental_strategy == "WGS") %>% droplevels()
tcgaGenusKrakenQCFiltWIS_WGS <- tcgaGenusKrakenQCFiltWIS[rownames(metadataSamplesAllQC_HiSeq_WGS),]

# Save data for ML
save(# Subset raw count data
  tcgaGenusKrakenQCFiltWIS_WGS,
  # Subset metadata
  metadataSamplesAllQC_HiSeq_WGS,
  file = "Interim_data/data_for_multiclass_ml_tcga_wgs24July22.RData")

# Scripts: SXXX

#----------------------------------------------------#
# Perform batch correction on WIS feature subset data
#----------------------------------------------------#

## VSNM
source("00-functions.R")
# Note: covDesignNorm had dimensions: 17625   211
batchCorrectedVSNMData <- vsnmFunctionTCGA(qcData = tcgaGenusKrakenQCFiltWIS,
                                           qcMetadata = metadataSamplesAllQC)
voomDataGenusKrakenQCFiltWIS <- batchCorrectedVSNMData$vdge_dataE
vsnmDataGenusKrakenQCFiltWIS <- batchCorrectedVSNMData$snmData

save(tcgaGenusKrakenQCFiltWIS,
     voomDataGenusKrakenQCFiltWIS,
     vsnmDataGenusKrakenQCFiltWIS,
     metadataSamplesAllQC,
     file = "Interim_data/data_vsnm_tcga_wis_features_subset_25July22.RData")

# NB: PCA plots and PVCA were run/generated on Barnacle under
# ~/projects/tcga/AA_Matters_Arising_Rebuttal

## Combatseq (using WGS-HiSeq data to correct for seqcenter efx)

#----------------------------------------------------#
# Save data for ML
#----------------------------------------------------#

save(# Subset raw count data
  tcgaGenusKrakenQCFiltWIS_HMS,
  tcgaGenusKrakenQCFiltWIS_BCM,
  tcgaGenusKrakenQCFiltWIS_MDA,
  tcgaGenusKrakenQCFiltWIS_WashU,
  tcgaGenusKrakenQCFiltWIS_Broad_WGS,
  tcgaGenusKrakenQCFiltWIS_UNC,
  tcgaGenusKrakenQCFiltWIS_CMS,
  tcgaGenusKrakenQCFiltWIS_Broad_RNA,
  # VSNM batch corrected data
  vsnmDataGenusKrakenQCFiltWIS,
  
  # Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  metadataSamplesAllQC_HiSeq_Broad_RNA,
  # Full metadata
  metadataSamplesAllQC,
  file = "Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy25July22.RData")

# Scripts: SXXX

#----------------------------------------------------#
# Plot PVCA results (originally run on Barnacle)
#----------------------------------------------------#
require(tibble)

load("Interim_data/pvcaVbRawNoVoomNoSNM_ExtendedFiltered.RData")
load("Interim_data/pvcaVoomNoSNM_ExtendedFiltered.RData")
load("Interim_data/pvcaSampleWithExpStrategySNM_ExtendedFiltered.RData")

pvcaData <- as.data.frame(rbind(pvcaVbRawNoVoomNoSNM_ExtendedFiltered_FA,
                                pvcaVoomNoSNM_ExtendedFiltered_FA,
                                pvcaSampleWithExpStrategySNM_ExtendedFiltered_FA)) %>%
  rownames_to_column("group")
pvcaData$group <- c("Raw Counts", "Voom", "Voom-SNM")
pvcaData.melted <- melt(pvcaData, id.vars = "group")
pvcaData.melted$group <- factor(pvcaData.melted$group, 
                                levels = c("Raw Counts","Voom","Voom-SNM"))

pvcaPlot <- ggplot(pvcaData.melted, aes(x = variable, y = value, fill = group)) + 
  geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
  geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(x = "Technical & Biological Effects",
       y = "Weighted average proportion variance",
       title = "PVCA of batch effect correction procedures\n(Only using WIS-overlapping genera)") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.3, 0.75)) + 
  scale_x_discrete(labels=c("sample_type" = "Sample Type", 
                            "disease_type" = "Disease Type", 
                            "data_submitting_center_label" = "Sequencing Center",
                            "platform" = "Sequencing Platform",
                            "experimental_strategy" = "Experimental strategy",
                            "tissue_source_site_label" = "Tissue Source Site",
                            "portion_is_ffpe" = "FFPE Fixation",
                            "resid" = "Residual \n(not explained by technical variation)")) +
  scale_fill_nejm(name = "Data types", labels = c("Raw Count Data", 
                                                  "Voom Normalized Data",
                                                  "Voom Normalized & SNM Corrected Data"))
ggsave(plot = pvcaPlot,
       filename = "Figures/PVCA-batch-correction-WIS-feature-subset.svg", 
       dpi = "retina",
       width = 18, height = 6, units = "in")

#----------------------------------------------------------#
# Examine % WIS-overlap with each version of decontamination
# Output: % of overlapping WIS features found in each decontam version
#----------------------------------------------------------#

load("../decontam4Datasets_Nov112019.RData", verbose = TRUE)

percentOverlapWIS <- function(df, wisFeat = wisGeneraBactKnownUniqueRev){
  krakenGenera <- gsub("k__.+\\.g__","",colnames(df))
  # Check overlap number
  print(sum(krakenGenera %in% wisFeat))
  print(sum(wisFeat %in% krakenGenera))
  
  # Intersect features
  krakenGeneraWISInteresected <- intersect(krakenGenera, wisFeat)
  
  # Print % of WIS-overlap
  print(sprintf("Percent overlap: %.2f", 100*length(krakenGeneraWISInteresected)/length(wisFeat)))
  
  # Create new data frames of intersected dataframes
  tcgaGenusKrakenQC <- df
  colnames(tcgaGenusKrakenQC) <- gsub("k__.+\\.g__","",colnames(tcgaGenusKrakenQC))
  tcgaGenusKrakenQCFiltWIS <- tcgaGenusKrakenQC[, krakenGeneraWISInteresected]
  
  return(tcgaGenusKrakenQCFiltWIS)
}

# NOTE: 218 total WIS genera in wisGeneraBactKnownUniqueRev
perOverlapWIS_full <- percentOverlapWIS(countsVbFinalNonzeroQC) # 184 / 84.40%
perOverlapWIS_likelyContamRemoved <- percentOverlapWIS(vbContaminatedDataQCLikelyContamRemoved) # 159 / 72.94%
perOverlapWIS_plateCenterContamRemoved <- percentOverlapWIS(vbContaminatedDataQCPlateCenterContamRemoved) # 129 / 59.17%
perOverlapWIS_allPutativeContamRemoved <- percentOverlapWIS(vbContaminatedDataQCAllPutativeContamRemoved) # 114 / 52.29%
perOverlapWIS_mostStringentContamRemoved <- percentOverlapWIS(vbContaminatedDataQCMostStringentContamRemoved) # 10 / 4.59%



