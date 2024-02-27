# 02-NR-Human-HPRC-reads.R
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
plotPrefix <- "K1C-HPRC"
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
countsVb_BIOM <- read_biom(biom_file = "Input_data/against-tcga-kraken1-processed/pangenome-adapter-filter-pese/genus_dbtcga_without_580_and_29_bac_and_vir_NONE_merged.biom")
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

#--------------Find overlap of WIS and Bins--------------#
binWISOverlap <- intersect(wisGeneraBactKnownUniqueRev,
                           binTumorGeneraKnownUnique)
length(binWISOverlap) # 69 (31.7% of 218)

countsVbFinalNonzeroBinsWIS <- tcgaGenusKrakenAll[, intersect(krakenGeneraWISInteresected,
                                                              binTumorGeneraKnownUnique)]
dim(countsVbFinalNonzeroBinsWIS) # 15512   69

#--------------Find overlap of BIO and Bins--------------#
load("Input_data/taxa_filtering_pipeline_13Oct23.RData",verbose = TRUE)

taxRS210_ff_combTaxaSpecies_UniqueGenera <- unique(taxRS210_ff_combTaxaSpecies$genus)
length(taxRS210_ff_combTaxaSpecies_UniqueGenera) # 442

binBIOOverlap <- intersect(taxRS210_ff_combTaxaSpecies_UniqueGenera,
                           binTumorGeneraKnownUnique)
length(binBIOOverlap) # 93 (31.7% of 218)

countsVbFinalNonzeroBinsBIO <- tcgaGenusKrakenAll[, intersect(colnames(tcgaGenusKrakenAll),
                                                              binBIOOverlap)]
dim(countsVbFinalNonzeroBinsBIO) # 15512   85

#--------------Stats for overlap of WIS and Bins--------------#
binTumorTaxaAndStatsKnownWIS <- binTumorTaxaAndStatsKnown %>%
  filter(genus %in% binWISOverlap)
dim(binTumorTaxaAndStatsKnownWIS) # 666 20

binTumorTaxaAndStatsKnownWIS_Complete <- binTumorTaxaAndStatsKnownWIS %>%
  filter(complete.cases(.)) %>% arrange(desc(size))
head(binTumorTaxaAndStatsKnownWIS_Complete)

binTumorTaxaAndStatsKnownWIS_Complete[,c("species","completeness","contamination","N50","size")] %>%
  write.csv("Outputs/binTumorTaxaAndStatsKnownWIS_Complete.csv",
            row.names = FALSE)
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
countsVbFinalNonzeroQCBinsBIO <- countsVbFinalNonzeroBinsBIO[rownames(metaDataFinalNonzeroQC),]

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

## Not currently converging
# ## VSNM -- BinsWIS feature subset data
# source("00-functions.R")
# # Note: covDesignNorm had dimensions: 15180    18
# batchCorrectedVSNMDataBinsWIS <- vsnmFunctionTCGA(qcData = countsVbFinalNonzeroQCBinsWIS,
#                                                qcMetadata = metaDataFinalNonzeroQC)
# voomDataGenusKrakenQCFiltBinsWIS <- batchCorrectedVSNMDataBinsWIS$vdge_dataE
# vsnmDataGenusKrakenQCFiltBinsWIS <- batchCorrectedVSNMDataBinsWIS$snmData

# Not currently converging
## VSNM -- BinsBIO feature subset data
source("00-functions.R")
batchCorrectedVSNMDataBinsBIO <- vsnmFunctionTCGA(qcData = countsVbFinalNonzeroQCBinsBIO,
                                                  qcMetadata = metaDataFinalNonzeroQC)
voomDataGenusKrakenQCFiltBinsBIO <- batchCorrectedVSNMDataBinsBIO$vdge_dataE
vsnmDataGenusKrakenQCFiltBinsBIO <- batchCorrectedVSNMDataBinsBIO$snmData

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
# Plot PVCA results (originally run on Barnacle)
#----------------------------------------------------#
require(tibble)
require(ggpubr)
require(ggsci)

load(paste0("Interim_data/",plotPrefix,"/pvcaData_Full.RData"), verbose = TRUE)
load(paste0("Interim_data/",plotPrefix,"/pvcaData_WIS.RData"), verbose = TRUE)

#--------------Full--------------#
pvcaDataFull <- as.data.frame(rbind(pvcaRaw_Full,
                                    pvcaVoom_Full,
                                    pvcaVSNM_Full)) %>%
  rownames_to_column("group")
pvcaDataFull$group <- c("Raw Counts", "Voom", "Voom-SNM")
pvcaDataFull.melted <- melt(pvcaDataFull, id.vars = "group")
pvcaDataFull.melted$group <- factor(pvcaDataFull.melted$group, 
                                    levels = c("Raw Counts","Voom","Voom-SNM"))

pvcaPlotFull <- ggplot(pvcaDataFull.melted, aes(x = variable, y = value, fill = group)) + 
  geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
  geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(x = "Technical & Biological Effects",
       y = "Weighted average proportion variance",
       title = "PVCA of batch effect correction procedures\n(Only using WIS-overlapping genera)") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.3, 0.8)) + 
  scale_x_discrete(labels=c("sample_type" = "Sample Type", 
                            "disease_type" = "Disease Type", 
                            "data_submitting_center_label" = "Sequencing Center",
                            "cgc_platform" = "Sequencing Platform",
                            "experimental_strategy" = "Experimental strategy",
                            "tissue_source_site_label" = "Tissue Source Site",
                            "portion_is_ffpe" = "FFPE Fixation",
                            "resid" = "Residual \n(not explained by technical variation)")) +
  scale_fill_nejm(name = "Data types", labels = c("Raw Count Data", 
                                                  "Voom Normalized Data",
                                                  "Voom Normalized & SNM Corrected Data"))
ggsave(plot = pvcaPlotFull,
       filename = paste0("Figures/",plotPrefix,"/pvcaPlot_Full.jpeg"), 
       dpi = "retina",
       width = 18, height = 6, units = "in")

#--------------WIS--------------#
pvcaDataWIS <- as.data.frame(rbind(pvcaRaw_WIS,
                                   pvcaVoom_WIS,
                                   pvcaVSNM_WIS)) %>%
  rownames_to_column("group")
pvcaDataWIS$group <- c("Raw Counts", "Voom", "Voom-SNM")
pvcaDataWIS.melted <- melt(pvcaDataWIS, id.vars = "group")
pvcaDataWIS.melted$group <- factor(pvcaDataWIS.melted$group, 
                                   levels = c("Raw Counts","Voom","Voom-SNM"))

pvcaPlotWIS <- ggplot(pvcaDataWIS.melted, aes(x = variable, y = value, fill = group)) + 
  geom_bar(aes(fill = group), position = "dodge", stat = "identity") + 
  geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(x = "Technical & Biological Effects",
       y = "Weighted average proportion variance",
       title = "PVCA of batch effect correction procedures\n(Only using WIS-overlapping genera)") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.3, 0.8)) + 
  scale_x_discrete(labels=c("sample_type" = "Sample Type", 
                            "disease_type" = "Disease Type", 
                            "data_submitting_center_label" = "Sequencing Center",
                            "cgc_platform" = "Sequencing Platform",
                            "experimental_strategy" = "Experimental strategy",
                            "tissue_source_site_label" = "Tissue Source Site",
                            "portion_is_ffpe" = "FFPE Fixation",
                            "resid" = "Residual \n(not explained by technical variation)")) +
  scale_fill_nejm(name = "Data types", labels = c("Raw Count Data", 
                                                  "Voom Normalized Data",
                                                  "Voom Normalized & SNM Corrected Data"))
ggsave(plot = pvcaPlotWIS,
       filename = paste0("Figures/",plotPrefix,"/pvcaPlot_WIS.jpeg"), 
       dpi = "retina",
       width = 18, height = 6, units = "in")

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

dim(countsVbFinalNonzeroQCForm) # 15180  1873
dim(vbDataBarnDFReconciledQC_Matched) # 15180  1873

#--------------------Corr--------------------#
require(ggpubr)
## By sample counts
cor.test(log10(rowSums(vbDataBarnDFReconciledQC_Matched)),
         log10(rowSums(countsVbFinalNonzeroQCForm)),
         method = "spearman") # 0.6820068

# Sanity check
all(rownames(countsVbFinalNonzeroQCForm) == metaDataFinalNonzeroQC$knightlabID) # TRUE
dfSampleCounts <- data.frame(Original = log10(rowSums(vbDataBarnDFReconciledQC_Matched)),
                             HPRC = log10(rowSums(countsVbFinalNonzeroQCForm)),
                             experimental_strategy = metaDataFinalNonzeroQC$experimental_strategy,
                             seqCenter = metaDataFinalNonzeroQC$data_submitting_center_label)

dfSampleCounts %>%
  ggscatter(x = "Original",
            y = "HPRC",
            xlab = "Original TCGA sample counts (log10)",
            ylab = "HPRC-depleted TCGA sample counts (log10)",
            color = "seqCenter",
            palette = "nejm",
            # legend = "right",
            # add = "reg.line", conf.int = TRUE,
            # add = "reg.line",
            # add.params = list(color = "blue"),
            # conf.int = TRUE,
            alpha = 0.07,
            # facet.by = "experimental_strategy",
            cor.coef = FALSE,
            cor.method = "spearman") +
  theme(aspect.ratio=1) +
  theme(text=element_text(size=14)) +
  stat_cor(method = "spearman",
           aes(color = seqCenter),
           show.legend = FALSE, label.y = seq(6,9,length.out = 8),
           cor.coef.name = "rho")
ggsave(filename = paste0("Figures/",plotPrefix,"/corr_by_sample_colored.jpeg"),
       width = 5.5, height = 5.5, units = "in", dpi = "retina")

## By microbe
# Remove zero-sum microbes
cor.test(log10(colSums(vbDataBarnDFReconciledQC_Matched)),
         log10(colSums(countsVbFinalNonzeroQCForm)),
         method = "spearman") # 0.7435244

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
            palette = c("black", "red"),
            xlab = "Original TCGA microbial counts (log10)",
            ylab = "HPRC-depleted TCGA microbial counts (log10)",
            # add = "reg.line",
            # add.params = list(color = "blue"),
            # conf.int = TRUE,
            alpha = 0.2,
            # facet.by = "experimental_strategy",
            cor.coef = FALSE) +
  theme(aspect.ratio=1) +
  theme(text=element_text(size=14)) +
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

wilcox.test(dfMicrobeCounts$Original[dfMicrobeCounts$WIS],
            dfMicrobeCounts$Original[!dfMicrobeCounts$WIS])

wilcox.test(dfMicrobeCounts$HPRC[dfMicrobeCounts$WIS],
            dfMicrobeCounts$HPRC[!dfMicrobeCounts$WIS])

#--------------------Compare genera--------------------#

# Number of genera in old data but not in new data --> 120
# Number of genera in new data but not in old data --> 1 (k__Viruses.o__Picornavirales.f__Picornaviridae.g__Teschovirus)
dim(vbDataBarnDFReconciledQC) # 17625  1993
dim(countsVbFinalNonzeroQC) # 15180  1874
sum(colnames(vbDataBarnDFReconciledQC) %in% colnames(countsVbFinalNonzeroQC)) # 1873
sum(colnames(countsVbFinalNonzeroQC) %in% colnames(vbDataBarnDFReconciledQC)) # 1873
colnames(countsVbFinalNonzeroQC)[which(!(colnames(countsVbFinalNonzeroQC) %in% colnames(vbDataBarnDFReconciledQC)))]
#----------------------------------------------------#
# Plot ML results -- heatmap
# Full vs WIS: PT, BDN, STN
#----------------------------------------------------#
require(tidyr)
require(pheatmap)
require(RColorBrewer)

abbreviationsTCGA_Allcancer <- read.csv("Input_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

csvFiles <- list.files(path = "Supporting_scripts/S03-ML-K1C-HPRC/", pattern="\\.csv$")
csvFilesFilt <- csvFiles[!grepl("_ALL_",csvFiles)]
csvList <- lapply(paste0("Supporting_scripts/S03-ML-K1C-HPRC/",csvFilesFilt), read.csv)
# names(csvList) <- csvFilesFilt
mlPerf_K1C_HPRC <- do.call(rbind, csvList) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=aucroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5)

# Format variable names
mlPerf_K1C_HPRC$datasetName[mlPerf_K1C_HPRC$datasetName == "vsnmDataGenusKrakenQCFilt"] <- "Full"
mlPerf_K1C_HPRC$datasetName[mlPerf_K1C_HPRC$datasetName == "vsnmDataGenusKrakenQCFiltWIS"] <- "WIS"
# mlPerf_K1C_HPRC$datasetName <- factor(mlPerf_K1C_HPRC$datasetName, levels = c("Full", "WIS"))
table(mlPerf_K1C_HPRC$datasetName)
table(mlPerf_K1C_HPRC$metadataName)

mlPerf_K1C_HPRC_VSNM <- mlPerf_K1C_HPRC %>%
  filter(metadataName == "metaDataFinalNonzeroQC") %>% droplevels()

## Heatmaps
heatmapPlot <- function(data=mlPerf_K1C_HPRC_VSNM, st, dt,meta=NULL,fontSize=16){
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

heatmapPlot2 <- function(data=mlPerf_K1C_HPRC_VSNM, st, dt,meta=NULL,fontSize=16){
  if(is.null(meta)){
    tmp <- data %>%
      filter(sampleType == st,
             datasetName == dt) %>%
      select(abbrev, AUROC, AUPR) %>%
      distinct(abbrev, .keep_all=TRUE) %>%
      arrange(abbrev) %>%
      right_join(abbreviationsTCGA_Allcancer, by = "abbrev") %>%
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
      right_join(abbreviationsTCGA_Allcancer, by = "abbrev") %>%
      arrange(abbrev) %>%
      column_to_rownames("abbrev")
  }
  breaksList <- seq(0, 1, length.out = 100)
  tmp %>% t() %>%
    pheatmap(cluster_rows = FALSE,
             cluster_cols = FALSE,
             border_color = "black",
             # color = colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(length(breaksList)),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
             breaks = breaksList,
             display_numbers = TRUE,
             number_format = "%.2f",
             na_col = "white",
             filename = paste0("Figures/",plotPrefix,"/heatmap_",
                               dt,"_",st,"__",paste(meta,sep="_"),".jpeg"),
             fontsize = fontSize,
             number_color = "black",
             width = 20,
             height = 1.5)
  tmp %>% summarise(meanAUROC = mean(AUROC, na.rm=TRUE),
                    meanAUPR = mean(AUPR,na.rm=TRUE)) %>% print()
}

#----------------Full----------------#
## PT
heatmapPlot2(st="Primary Tumor", dt="Full", meta="metaDataFinalNonzeroQC")
## BDN
heatmapPlot2(st="Blood Derived Normal", dt="Full", meta="metaDataFinalNonzeroQC")
## STN
heatmapPlot2(st="Primary Tumor vs Solid Tissue Normal", dt="Full", meta="metaDataFinalNonzeroQC")

#----------------WIS----------------#
## PT
heatmapPlot2(st="Primary Tumor", dt="WIS", meta="metaDataFinalNonzeroQC")
## BDN
heatmapPlot2(st="Blood Derived Normal", dt="WIS", meta="metaDataFinalNonzeroQC")
## STN
heatmapPlot2(st="Primary Tumor vs Solid Tissue Normal", dt="WIS", meta="metaDataFinalNonzeroQC")

## Combined
breaksList <- seq(0, 1, length.out = 100)
mlPerf_K1C_HPRC_VSNM %>%
  filter(metadataName == "metaDataFinalNonzeroQC",
         sampleType == "Blood Derived Normal") %>%
  select(abbrev, datasetName, AUROC, AUPR) %>%
  arrange(abbrev) %>%
  pivot_wider(names_from = "datasetName", id_cols = "abbrev", values_from = c("AUROC","AUPR")) %>%
  column_to_rownames("abbrev") %>%
  relocate(AUPR_Full, .after = AUROC_Full) %>%
  pheatmap(cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = "black",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           breaks = breaksList,
           display_numbers = TRUE,
           number_format = "%.2f",
           na_col = "white",
           filename = paste0("Figures/",plotPrefix,"/heatmap_BDN_All.jpeg"),
           fontsize = 14,
           number_color = "black",
           width = 4,
           height = 8)

#----------------BinsWIS----------------#

mlPerf_K1C_HPRC_BinsWIS <- mlPerf_K1C_HPRC %>%
  filter(grepl("BinsWIS",datasetName)) %>%
  mutate(metadataName = gsub("metaDataFinalNonzeroQC_HiSeq_","",metadataName)) %>%
  mutate(seqCenter = metadataName)

mlPerf_K1C_HPRC_BinsWIS %>%
  mutate(sampleType = gsub("Primary Tumor vs Solid Tissue Normal","Tumor vs Normal",sampleType)) %>%
  select(AUROC, AUPR, abbrev,diseaseType,
         sampleType,datasetName, seqCenter) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor","Blood Derived Normal","Tumor vs Normal"))) %>%
  ggboxplot(x = "seqCenter",
            y = "AUROC",
            fill = "datasetName",
            legend = "none",
            # add = "jitter",
            # add.params = list(alpha=0.3),
            order = c("HMS","MDA","WashU","Broad_WGS","BCM","CMS","UNC"),
            palette = "nejm") +
  facet_grid( ~ sampleType, scales = "free") +
  rotate_x_text(30) +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1), limits = c(0.4,1))
ggsave(filename = paste0("Figures/",plotPrefix,"/mlPerf_BinsWIS_Boxplot_AUROC_23Feb24.jpeg"),
       dpi = "retina", units = "in", width = 8, height = 5)

mlPerf_K1C_HPRC_BinsWIS %>%
  mutate(sampleType = gsub("Primary Tumor vs Solid Tissue Normal","Tumor vs Normal",sampleType)) %>%
  select(AUROC, AUPR, abbrev,diseaseType,
         sampleType,datasetName, seqCenter) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor","Blood Derived Normal","Tumor vs Normal"))) %>%
  ggboxplot(x = "seqCenter",
            y = "AUPR",
            fill = "datasetName",
            legend = "none",
            # add = "jitter",
            # add.params = list(alpha=0.3),
            order = c("HMS","MDA","WashU","Broad_WGS","BCM","CMS","UNC"),
            palette = "nejm") +
  facet_grid( ~ sampleType, scales = "free") +
  rotate_x_text(30) +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1))
ggsave(filename = paste0("Figures/",plotPrefix,"/mlPerf_BinsWIS_Boxplot_AUPR_23Feb24.jpeg"),
       dpi = "retina", units = "in", width = 8, height = 5)

mlPerf_K1C_HPRC_BinsWIS %>%
  mutate(sampleType = gsub("Primary Tumor vs Solid Tissue Normal","Tumor vs Normal",sampleType)) %>%
  select(AUROC, AUPR, abbrev,diseaseType,
         sampleType,datasetName, seqCenter) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor","Blood Derived Normal","Tumor vs Normal"))) %>%
  group_by(sampleType) %>%
  # group_by(seqCenter, sampleType) %>%
  summarise(meanAUROC = mean(AUROC),
            meanAUPR = mean(AUPR)) %>%
  arrange(sampleType)

#----------------------------------------------------#
# Plot ML results -- heatmap
# Full vs WIS: Stage I vs IV
#----------------------------------------------------#

require(tidyr)
require(pheatmap)
abbreviationsTCGA_Allcancer <- read.csv("Input_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

## Load data
mlPerf_K1C_HPRC_StageIvsIV <- read.csv("Supporting_scripts/S06-ML-K1C-HPRC-Stage/perfML_Human_HPRC_StageIvsIV_ALL_22Feb24.csv", 
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
mlPerf_K1C_HPRC_StageIvsIV$datasetName[mlPerf_K1C_HPRC_StageIvsIV$datasetName == "vsnmDataGenusKrakenQCFilt_Path"] <- "Full"
mlPerf_K1C_HPRC_StageIvsIV$datasetName[mlPerf_K1C_HPRC_StageIvsIV$datasetName == "vsnmDataGenusKrakenQCFiltWIS_Path"] <- "WIS"
# mlPerf_K1C_HPRC_StageIvsIV$datasetName <- factor(mlPerf_K1C_HPRC_StageIvsIV$datasetName, levels = c("Full", "WIS"))
table(mlPerf_K1C_HPRC_StageIvsIV$datasetName)
table(mlPerf_K1C_HPRC_StageIvsIV$metadataName)

## Full
heatmapPlot2(data = mlPerf_K1C_HPRC_StageIvsIV, st="Stage I vs IV", 
             dt="Full", meta="metaDataFinalNonzeroQCPath")
## WIS
heatmapPlot2(data = mlPerf_K1C_HPRC_StageIvsIV, st="Stage I vs IV", 
             dt="WIS", meta="metaDataFinalNonzeroQCPath")

#----------------------------------------------------#
# Plot ML results -- heatmap
# Full vs WIS: BDN early stage, BDN noMut
#----------------------------------------------------#

require(tidyr)
require(pheatmap)
require(RColorBrewer)

breaksList <- seq(0, 1, length.out = 100)

abbreviationsTCGA_Allcancer <- read.csv("Input_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

## Load data
mlPerf_K1C_HPRC_BDNStage <- read.csv("Supporting_scripts/S06-ML-K1C-HPRC-Stage/perfML_Human_HPRC_BDN_EarlyStage_noMut_ALL_22Feb24.csv", 
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
mlPerf_K1C_HPRC_BDNStage$datasetName <- ifelse(grepl("WIS_",mlPerf_K1C_HPRC_BDNStage$datasetName),
                                               "WIS","Full")
table(mlPerf_K1C_HPRC_BDNStage$datasetName)
table(mlPerf_K1C_HPRC_BDNStage$metadataName)

## Heatmaps
heatmapPlot3 <- function(data=mlPerf_K1C_HPRC_VSNM, st, dt,meta=NULL,fontSize=16){
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
  tmp %>%
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

#----------------BDN early stage----------------#
mlPerf_K1C_HPRC_BDNStage %>%
  filter(metadataName == "metaDataFinalNonzeroQCPath_BDN_EarlyStage") %>%
  select(abbrev, datasetName, AUROC, AUPR) %>%
  arrange(abbrev) %>%
  pivot_wider(names_from = "datasetName", values_from = c("AUROC","AUPR")) %>%
  column_to_rownames("abbrev") %>%
  relocate(AUPR_Full, .after = AUROC_Full) %>%
  pheatmap(cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = "black",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           breaks = breaksList,
           display_numbers = TRUE,
           number_format = "%.2f",
           na_col = "white",
           filename = paste0("Figures/",plotPrefix,"/heatmap_BDN_EarlyStage.jpeg"),
           fontsize = 14,
           number_color = "black",
           width = 4,
           height = 4)

mlPerf_K1C_HPRC_BDNStage %>%
  filter(metadataName == "metaDataFinalNonzeroQCPath_BDN_EarlyStage") %>%
  select(abbrev, datasetName, AUROC, AUPR) %>%
  arrange(abbrev) %>%
  group_by(datasetName) %>% 
  summarise(meanAUROC = mean(AUROC))

#----------------BDN noMut----------------#
mlPerf_K1C_HPRC_BDNStage %>%
  filter(metadataName == "metaDataFinalNonzeroQC_noMutGuardant") %>%
  select(abbrev, datasetName, AUROC, AUPR) %>%
  arrange(abbrev) %>%
  pivot_wider(names_from = "datasetName", values_from = c("AUROC","AUPR")) %>%
  column_to_rownames("abbrev") %>%
  relocate(AUPR_Full, .after = AUROC_Full) %>%
  pheatmap(cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = "black",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           breaks = breaksList,
           display_numbers = TRUE,
           number_format = "%.2f",
           na_col = "white",
           filename = paste0("Figures/",plotPrefix,"/heatmap_BDN_noMutGuardant.jpeg"),
           fontsize = 14,
           number_color = "black",
           width = 4,
           height = 4)

mlPerf_K1C_HPRC_BDNStage %>%
  filter(metadataName == "metaDataFinalNonzeroQC_noMutGuardant") %>%
  select(abbrev, datasetName, AUROC, AUPR) %>%
  arrange(abbrev) %>%
  group_by(datasetName) %>% 
  summarise(meanAUROC = mean(AUROC))

mlPerf_K1C_HPRC_BDNStage %>%
  filter(metadataName == "metaDataFinalNonzeroQC_noMutFoundation") %>%
  select(abbrev, datasetName, AUROC, AUPR) %>%
  arrange(abbrev) %>%
  pivot_wider(names_from = "datasetName", values_from = c("AUROC","AUPR")) %>%
  column_to_rownames("abbrev") %>%
  relocate(AUPR_Full, .after = AUROC_Full) %>%
  pheatmap(cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = "black",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           breaks = breaksList,
           display_numbers = TRUE,
           number_format = "%.2f",
           na_col = "white",
           filename = paste0("Figures/",plotPrefix,"/heatmap_BDN_noMutFoundation.jpeg"),
           fontsize = 14,
           number_color = "black",
           width = 4,
           height = 2.5)

mlPerf_K1C_HPRC_BDNStage %>%
  filter(metadataName == "metaDataFinalNonzeroQC_noMutFoundation") %>%
  select(abbrev, datasetName, AUROC, AUPR) %>%
  arrange(abbrev) %>%
  group_by(datasetName) %>% 
  summarise(meanAUROC = mean(AUROC))

## Early stage
heatmapPlot(data = mlPerf_K1C_HPRC_BDNStage, st = "Blood Derived Normal", 
            dt = "Full", meta = "metaDataFinalNonzeroQCPath_BDN_EarlyStage")
## noMut Guardant
heatmapPlot2(data = mlPerf_K1C_HPRC_BDNStage, st = "Blood Derived Normal", 
             dt = "Full", meta = "metaDataFinalNonzeroQC_noMutGuardant")
## noMut Foundation
heatmapPlot2(data = mlPerf_K1C_HPRC_BDNStage, st = "Blood Derived Normal", 
             dt = "Full", meta = "metaDataFinalNonzeroQC_noMutFoundation")

#----------------WIS: Early stage and noMut BDN----------------#
## Early stage
heatmapPlot2(data = mlPerf_K1C_HPRC_BDNStage, st = "Blood Derived Normal", 
             dt = "WIS", meta = "metaDataFinalNonzeroQCPath_BDN_EarlyStage")
## noMut Guardant
heatmapPlot2(data = mlPerf_K1C_HPRC_BDNStage, st = "Blood Derived Normal", 
             dt = "WIS", meta = "metaDataFinalNonzeroQC_noMutGuardant")
## noMut Foundation
heatmapPlot2(data = mlPerf_K1C_HPRC_BDNStage, st = "Blood Derived Normal", 
             dt = "WIS", meta = "metaDataFinalNonzeroQC_noMutFoundation")
#----------------------------------------------------#
# Plot ML results: WGS vs RNA
#----------------------------------------------------#
require(tidyr)
require(pheatmap)

abbreviationsTCGA_Allcancer <- read.csv("Input_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

mlPerf_K1C_HPRC_WGSvsRNA <- read.csv("Supporting_scripts/S07-ML-K1C-HPRC-WGSvsRNA/perfML_Human_HPRC_WGSvsRNA_ALL_22Feb24.csv", 
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
mlPerf_K1C_HPRC_WGSvsRNA$datasetName[mlPerf_K1C_HPRC_WGSvsRNA$datasetName == "vsnmDataGenusKrakenQCFilt_WGS"] <- "Full WGS"
mlPerf_K1C_HPRC_WGSvsRNA$datasetName[mlPerf_K1C_HPRC_WGSvsRNA$datasetName == "vsnmDataGenusKrakenQCFilt_RNA"] <- "Full RNA"
mlPerf_K1C_HPRC_WGSvsRNA$datasetName[mlPerf_K1C_HPRC_WGSvsRNA$datasetName == "vsnmDataGenusKrakenQCFiltWIS_WGS"] <- "WIS WGS"
mlPerf_K1C_HPRC_WGSvsRNA$datasetName[mlPerf_K1C_HPRC_WGSvsRNA$datasetName == "vsnmDataGenusKrakenQCFiltWIS_RNA"] <- "WIS RNA"
mlPerf_K1C_HPRC_WGSvsRNA$metadataName[mlPerf_K1C_HPRC_WGSvsRNA$metadataName == "metaDataFinalNonzeroQC_WGS"] <- "WGS"
mlPerf_K1C_HPRC_WGSvsRNA$metadataName[mlPerf_K1C_HPRC_WGSvsRNA$metadataName == "metaDataFinalNonzeroQC_RNA"] <- "RNA"
mlPerf_K1C_HPRC_WGSvsRNA$datasetName <- factor(mlPerf_K1C_HPRC_WGSvsRNA$datasetName,
                                               levels = c("Full WGS",
                                                          "Full RNA",
                                                          "WIS WGS",
                                                          "WIS RNA"))
table(mlPerf_K1C_HPRC_WGSvsRNA$datasetName)
table(mlPerf_K1C_HPRC_WGSvsRNA$metadataName)

#----------------Full----------------#
## PT
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor", 
             dt="Full WGS")
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor", 
             dt="Full RNA")
## STN
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor vs Solid Tissue Normal", 
             dt="Full WGS")
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor vs Solid Tissue Normal", 
             dt="Full RNA")

#----------------WIS----------------#
## PT
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor", 
             dt="WIS WGS")
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor", 
             dt="WIS RNA")
## STN
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor vs Solid Tissue Normal", 
             dt="WIS WGS")
heatmapPlot2(mlPerf_K1C_HPRC_WGSvsRNA, st="Primary Tumor vs Solid Tissue Normal", 
             dt="WIS RNA")

#----------------------------------------------------#
# Plot ML results -- actual vs negative control
#----------------------------------------------------#
require(rstatix)
# Use mlPerf_K1C_HPRC object (computed above)

mlPerf_K1C_HPRC_Shu <- read.csv("Supporting_scripts/S08-ML-K1C-HPRC-shuffled-counts/perfML_Human_HPRC_Full_WIS_Bins_Shuffled_ALL_22Feb24.csv", 
                                stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=aucroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5) %>%
  mutate(datasetName = gsub("vsnmDataGenusKrakenQCFiltWIS","VSNM_WIS",datasetName)) %>%
  mutate(datasetName = gsub("vsnmDataGenusKrakenQCFilt","VSNM_Full",datasetName))

mlPerf_K1C_HPRC_Scr <- read.csv("Supporting_scripts/S09-ML-K1C-HPRC-scrambled-metadata/perfML_Human_HPRC_Full_WIS_Bins_scrambled_ALL_22Feb24.csv", 
                                stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=aucroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5) %>%
  mutate(datasetName = gsub("vsnmDataGenusKrakenQCFiltWIS","VSNM_WIS",datasetName)) %>%
  mutate(datasetName = gsub("vsnmDataGenusKrakenQCFilt","VSNM_Full",datasetName))

mlPerf_K1C_HPRC_Scr_Shu <- mlPerf_K1C_HPRC %>%
  mutate(datasetName = gsub("^WIS","VSNM_WIS",datasetName)) %>%
  mutate(datasetName = gsub("^Full","VSNM_Full",datasetName)) %>%
  bind_rows(mlPerf_K1C_HPRC_Shu,
            mlPerf_K1C_HPRC_Scr) %>%
  mutate(datasetName = gsub("countsVbFinalNonzeroQCBins_HiSeq","Raw_Bins",datasetName)) %>%
  mutate(datasetName = gsub("countsVbFinalNonzeroQCBinsWIS_HiSeq","Raw_BinsWIS",datasetName)) %>%
  mutate(datasetName = gsub("countsVbFinalNonzeroQCWIS_HiSeq","Raw_WIS",datasetName)) %>%
  mutate(datasetName = gsub("countsVbFinalNonzeroQC_HiSeq","Raw_Full",datasetName)) %>%
  mutate(metadataName = gsub("metaDataFinalNonzeroQC_HiSeq_","",metadataName)) %>%
  mutate(metadataName = gsub("_HMS|_BCM|_MDA|_WashU|_UNC|_CMS","",metadataName)) %>%
  mutate(metadataName = gsub("_scrambled","_Scr",metadataName)) %>%
  mutate(metadataName = gsub("_shuffled","_Shu",metadataName)) %>%
  mutate(metadataName = case_when(
    (metadataName == "metaDataFinalNonzeroQC") & (datasetName=="VSNM_WIS") ~ "WIS",
    (metadataName == "metaDataFinalNonzeroQC_Scr") & (datasetName=="VSNM_WIS") ~ "WIS_Scr",
    (metadataName == "metaDataFinalNonzeroQC_Shu") & (datasetName=="VSNM_WIS") ~ "WIS_Shu",
    (metadataName == "metaDataFinalNonzeroQC") & (datasetName=="VSNM_Full") ~ "Full",
    (metadataName == "metaDataFinalNonzeroQC_Scr") & (datasetName=="VSNM_Full") ~ "Full_Scr",
    (metadataName == "metaDataFinalNonzeroQC_Shu") & (datasetName=="VSNM_Full") ~ "Full_Shu",
    .default = metadataName
  )) %>%
  mutate(sampleType = gsub("Primary Tumor vs Solid Tissue Normal","PT vs STN",sampleType)) %>%
  mutate(sampleType = gsub("Primary Tumor","PT",sampleType)) %>%
  mutate(sampleType = gsub("Blood Derived Normal","BDN",sampleType)) %>%
  select(abbrev,sampleType,AUROC,AUPR,datasetName,metadataName)

table(mlPerf_K1C_HPRC_Scr_Shu$datasetName)
table(mlPerf_K1C_HPRC_Scr_Shu$metadataName)

plotActualVsCtrls_VSNM <- function(dt = "VSNM_Full"){
  plotDataWilcoxAUROC <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    group_by(sampleType) %>%
    wilcox_test(AUROC ~ metadataName) %>%
    adjust_pvalue(method = "BH") %>%
    p_round()
  
  plotDataWilcoxAUPR <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    group_by(sampleType) %>%
    wilcox_test(AUPR ~ metadataName) %>%
    adjust_pvalue(method = "BH") %>%
    p_round()
  
  pAUROC <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    ggboxplot(x = "metadataName",
              y = "AUROC",
              facet.by = "sampleType",
              fill = "metadataName",
              ylim = c(0,1.3),
              legend = "none",
              palette = "nejm",
              notch = TRUE) +
    ggbeeswarm::geom_quasirandom(alpha=0.1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  
  pAUPR <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    ggboxplot(x = "metadataName",
              y = "AUPR",
              facet.by = "sampleType",
              fill = "metadataName",
              ylim = c(0,1.3),
              legend = "none",
              palette = "nejm",
              notch = TRUE) +
    ggbeeswarm::geom_quasirandom(alpha=0.1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  
  plotDataWilcoxAUROC <- plotDataWilcoxAUROC %>% add_xy_position(x="metadataName")
  plotDataWilcoxAUPR <- plotDataWilcoxAUPR %>% add_xy_position(x="metadataName")
  pAUROC + stat_pvalue_manual(plotDataWilcoxAUROC, y.position = c(1.1,1.2,1.3),
                              size = 3) + rotate_x_text(30)
  ggsave(filename = paste0("Figures/",plotPrefix,"/auroc_actual_vs_ctrl_",dt,".jpeg"),
         dpi = "retina", units = "in", width = 3.5, height = 4.5)
  pAUPR + stat_pvalue_manual(plotDataWilcoxAUPR,  y.position = c(1.1,1.2,1.3),
                             size=3) + rotate_x_text(30)
  ggsave(filename = paste0("Figures/",plotPrefix,"/aupr_actual_vs_ctrl_",dt,".jpeg"),
         dpi = "retina", units = "in", width = 3.5, height = 4.5)
}

plotActualVsCtrls_SeqCenter <- function(dt = "VSNM_Full"){
  plotDataWilcoxAUROC <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    # group_by(sampleType) %>%
    wilcox_test(AUROC ~ metadataName) %>%
    adjust_pvalue(method = "BH") %>%
    p_round()
  
  plotDataWilcoxAUPR <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    # group_by(sampleType) %>%
    wilcox_test(AUPR ~ metadataName) %>%
    adjust_pvalue(method = "BH") %>%
    p_round()
  
  pAUROC <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    ggboxplot(x = "metadataName",
              y = "AUROC",
              # facet.by = "sampleType",
              fill = "metadataName",
              ylim = c(0,1.3),
              legend = "none",
              palette = "nejm",
              notch = FALSE) +
    ggbeeswarm::geom_quasirandom(alpha=0.1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  
  pAUPR <- mlPerf_K1C_HPRC_Scr_Shu %>%
    filter(grepl(dt,datasetName)) %>%
    ggboxplot(x = "metadataName",
              y = "AUPR",
              # facet.by = "sampleType",
              fill = "metadataName",
              ylim = c(0,1.3),
              legend = "none",
              palette = "nejm",
              notch = FALSE) +
    ggbeeswarm::geom_quasirandom(alpha=0.1) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  
  plotDataWilcoxAUROC <- plotDataWilcoxAUROC %>% add_xy_position(x="metadataName")
  plotDataWilcoxAUPR <- plotDataWilcoxAUPR %>% add_xy_position(x="metadataName")
  pAUROC + stat_pvalue_manual(plotDataWilcoxAUROC, y.position = c(1.1,1.2,1.3),
                              size = 3) + rotate_x_text(30)
  ggsave(filename = paste0("Figures/",plotPrefix,"/auroc_actual_vs_ctrl_",dt,".jpeg"),
         dpi = "retina", units = "in", width = 3.5, height = 4.5)
  pAUPR + stat_pvalue_manual(plotDataWilcoxAUPR,  y.position = c(1.1,1.2,1.3),
                             size=3) + rotate_x_text(30)
  ggsave(filename = paste0("Figures/",plotPrefix,"/aupr_actual_vs_ctrl_",dt,".jpeg"),
         dpi = "retina", units = "in", width = 3.5, height = 4.5)
}

table(mlPerf_K1C_HPRC_Scr_Shu$datasetName)
# VSNM
plotActualVsCtrls_VSNM(dt = "VSNM_Full")
plotActualVsCtrls_VSNM(dt = "VSNM_WIS")

# # Full -- not run
# plotActualVsCtrls_SeqCenter(dt = "Raw_Full_HMS")
# plotActualVsCtrls_SeqCenter(dt = "Raw_Full_BCM")
# plotActualVsCtrls_SeqCenter(dt = "Raw_Full_MDA")
# plotActualVsCtrls_SeqCenter(dt = "Raw_Full_WashU")
# plotActualVsCtrls_SeqCenter(dt = "Raw_Full_UNC")
# plotActualVsCtrls_SeqCenter(dt = "Raw_Full_CMS")
# 
# # WIS -- not run
# plotActualVsCtrls_SeqCenter(dt = "Raw_WIS_HMS")
# plotActualVsCtrls_SeqCenter(dt = "Raw_WIS_BCM")
# plotActualVsCtrls_SeqCenter(dt = "Raw_WIS_MDA")
# plotActualVsCtrls_SeqCenter(dt = "Raw_WIS_WashU")
# plotActualVsCtrls_SeqCenter(dt = "Raw_WIS_UNC")
# plotActualVsCtrls_SeqCenter(dt = "Raw_WIS_CMS")

#----------------------------------------------------#
# Original ecological validation analyses
#----------------------------------------------------#
require(ggpubr)
require(ggsci)

load(paste0("Interim_data/",plotPrefix,"/data_vsnm_tcga_full_wis_bins_features_subset_20Feb24.RData"),
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
hpvPancancerData <- data.frame(HPV = vsnmDataGenusKrakenQCFilt[rownames(metaDataFinalNonzeroQC),
                                                               "k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus"])
hpvPancancerCombined <- droplevels(cbind(hpvPancancerMeta, hpvPancancerData))

interactVec <- as.character(interaction(hpvPancancerCombined$investigation, 
                                        hpvPancancerCombined$hpv_status, sep = " "))
interactVec[which(is.na(interactVec))] <- as.character(hpvPancancerCombined$investigation[which(is.na(interactVec))])
interactVec <- factor(interactVec)
hpvPancancerCombined$hpvInteract <- interactVec

# Trying out HNSC with this metadata
hpvHNSCComparisons <- list( c("TCGA-HNSC Positive", "TCGA-HNSC Negative"))
hpvPancancerCombined %>%
  filter((sample_type %in% c("Blood Derived Normal", "Primary Tumor")) &
           !(hpv_status %in% c("Indeterminate")) &
           (investigation == "TCGA-HNSC")) %>%
  ggboxplot(x = "hpvInteract", y = "HPV", 
            # color = "sample_type",
            add = "jitter",
            facet.by = "sample_type",
            # palette = pal_nejm(),
            xlab = "Clinical HPV Status", ylab = "SNM Normalized Abundance",
            ylim = c(7, 20),
            title = "Pancancer Comparison of Alphapapillomavirus Genus Abundance in HNSC",
            # legend = "right",
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels=c("TCGA-HNSC Negative" = "Negative",
                            "TCGA-HNSC Positive" = "Positive")) +
  scale_color_nejm() +
  rotate_x_text(angle = 30) +
  stat_summary(geom = "text", 
               fun.data = function(x){c(y = 7, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0.75)) +
  stat_compare_means(comparisons = hpvHNSCComparisons, label.y = 19,
                     label = "p.format") #-> p # Add pairwise comparisons p-value
ggsave(filename = paste0("Figures/",plotPrefix,"/ecologicalVal_HPV_HNSC_bigrquery.jpeg"),
       width = 4, height = 4.5, units = "in", dpi = "retina")

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
            # ylim = c(-1, 21),
            # order = c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal"),
            legend.title = "Sample Type") +
  rotate_x_text(angle = 30) +
  theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
  stat_summary(geom = "text", 
               fun.data = function(x){c(y = 7.5, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0.75)) +
  stat_compare_means(comparisons = coadFusoGIComp, 
                     label.y = c(20,20,20),
                     #size = 8,
                     label = "p.format") # -> fusoPlotGI # Add pairwise comparisons p-value
ggsave(filename = paste0("Figures/",plotPrefix,"/ecologicalVal_Fuso_PT_STN.jpeg"),
       width = 5, height = 6, units = "in", dpi = "retina")

coadFusoCombined %>%
  # filter(experimental_strategy %in% c("RNA-Seq")) %>%
  # filter(pathologic_stage_label %in% c("Stage IA","Stage IIA","Stage IIIA","Stage IVA")) %>%
  filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
  # filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  # filter(sample_type %in% c("Primary Tumor")) %>%
  filter(disease_type %in% c("Colon Adenocarcinoma", "Stomach Adenocarcinoma", "Esophageal Carcinoma", "Rectum Adenocarcinoma")) %>%
  # filter(disease_type %in% c("Rectum Adenocarcinoma")) %>%
  mutate(sample_type = gsub("Solid Tissue Normal","STN",sample_type)) %>%
  mutate(sample_type = gsub("Primary Tumor","PT",sample_type)) %>%
  mutate(sample_type = gsub("Blood Derived Normal","BDN",sample_type)) %>%
  ggboxplot(x = "sample_type", y = "Fuso", 
            # color = "GIBoolean",
            add = "jitter",
            add.params = list(alpha = 0.1),
            palette = "lancet",
            facet.by = "investigation",
            nrow = 1,
            xlab = "Cancer Type", ylab = "SNM Normalized Abundance", 
            title = "Comparison of Fusobacterium Genus Abundance",
            legend = "right",
            # ylim = c(-1, 21),
            order = c("STN", "PT", "BDN"),
            legend.title = "Sample Type") +
  # rotate_x_text(angle = 30) +
  theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
  stat_summary(geom = "text", 
               fun.data = function(x){c(y = 7.5, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0.75)) +
  stat_compare_means(comparisons = list(c("STN","PT"),
                                        c("PT","BDN"),
                                        c("STN","BDN")), 
                     label.y = c(20,22,24),
                     #size = 8,
                     label = "p.format") # -> fusoPlotGI # Add pairwise comparisons p-value
ggsave(filename = paste0("Figures/",plotPrefix,"/ecologicalVal_Fuso_4CTs_with_STs.jpeg"),
       width = 7, height = 5, units = "in", dpi = "retina")

#----------------------------------------------------------#
# Calculate alpha and beta diversity per batch using phyloseq
#----------------------------------------------------------#
require(ggsci)
require(ggpubr)

source("00-functions.R") # for alphaBetaFXN() and runAlphaBetaSeqCenter()
# WIS
runAlphaBetaSeqCenter(metaString = "metaDataFinalNonzeroQC_HiSeq",
                      countString = "countsVbFinalNonzeroQCWIS_HiSeq",
                      dataStringInput = "K1C_HPRC_WIS",
                      useTaxTableInput = FALSE)
ab_UNC_WIS <- alphaOnlyFXN(metaData = eval(as.name(paste0("metaDataFinalNonzeroQC_HiSeq","_UNC"))),
                           countData = eval(as.name(paste0("countsVbFinalNonzeroQCWIS_HiSeq","_UNC"))),
                           alphaDivType = "Observed",
                           useTaxTable = FALSE,
                           dataString = paste0("K1C_HPRC_WIS","_UNC"),
                           alphaPlotHeight = 2.5,
                           alphaPlotWidth = 8,
                           ptOnlyFlag = TRUE)

# BinsWIS
runAlphaBetaSeqCenter(metaString = "metaDataFinalNonzeroQC_HiSeq",
                      countString = "countsVbFinalNonzeroQCBinsWIS_HiSeq",
                      dataStringInput = "K1C_HPRC_BinsWIS",
                      useTaxTableInput = FALSE)

ab_UNC_BinsWIS <- alphaOnlyFXN(metaData = eval(as.name(paste0("metaDataFinalNonzeroQC_HiSeq","_UNC"))),
                               countData = eval(as.name(paste0("countsVbFinalNonzeroQCBinsWIS_HiSeq","_UNC"))),
                               alphaDivType = "Observed",
                               useTaxTable = FALSE,
                               dataString = paste0("K1C_HPRC_BinsWIS","_UNC"),
                               alphaPlotHeight = 2.5,
                               alphaPlotWidth = 8,
                               ptOnlyFlag = TRUE)
# # Full
# runAlphaBetaSeqCenter(metaString = "metaDataFinalNonzeroQC_HiSeq",
#                       countString = "countsVbFinalNonzeroQC_HiSeq",
#                       dataStringInput = "K1C_HPRC_Full",
#                       useTaxTableInput = FALSE)

#----------------------------------------------------#
# Export per-batch data to run Qiime2
#----------------------------------------------------#
require(biomformat)
require(rhdf5)

## Create and write taxa file
taxaFileWIS <- data.frame(`Feature ID` = colnames(countsVbFinalNonzeroQCWIS_HiSeq),
                          Taxon = colnames(countsVbFinalNonzeroQCWIS_HiSeq),
                          check.names = FALSE)
write.table(taxaFileWIS,
            file = "Qiime_data_and_scripts/Qiime_input_data/wis-taxa.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

source("00-functions.R") # for export2Qiime()
# WIS
export2Qiime(metaString = "metaDataFinalNonzeroQC_HiSeq",
             countString = "countsVbFinalNonzeroQCWIS_HiSeq",
             dataString = "K1C_HPRC_WIS")

# BinsWIS
export2Qiime(metaString = "metaDataFinalNonzeroQC_HiSeq",
             countString = "countsVbFinalNonzeroQCBinsWIS_HiSeq",
             dataString = "K1C_HPRC_BinsWIS")

# # Full
# export2Qiime(metaString = "metaDataFinalNonzeroQC_HiSeq",
#              countString = "countsVbFinalNonzeroQC_HiSeq",
#              dataString = "K1C_HPRC_Full")

#----------------------------------------------------#
# Calculate per-center differential abundances
#----------------------------------------------------#

source("00-functions.R") # for runAncomBC_1VsAll_OGUs()
# WIS
runAncomBC_1VsAll_OGUs(metaString = "metaDataFinalNonzeroQC_HiSeq",
                       countString = "countsVbFinalNonzeroQCWIS_HiSeq",
                       dataString = "K1C_HPRC_WIS",
                       taxTable = NA,
                       makeTaxTable = TRUE,
                       qvalCutoff = 0.05,
                       sampleTypes = c("Primary Tumor","Blood Derived Normal"),
                       SeqCenters = c("HMS","BCM","MDA","WashU","Broad_WGS","CMS","UNC"),
                       taxaPlotLabel = "genus")

# BinsWIS
runAncomBC_1VsAll_OGUs(metaString = "metaDataFinalNonzeroQC_HiSeq",
                       countString = "countsVbFinalNonzeroQCBinsWIS_HiSeq",
                       dataString = "K1C_HPRC_BinsWIS",
                       taxTable = NA,
                       makeTaxTable = TRUE,
                       qvalCutoff = 0.05,
                       sampleTypes = c("Primary Tumor","Blood Derived Normal"),
                       SeqCenters = c("HMS","BCM","MDA","WashU","Broad_WGS","CMS","UNC"),
                       taxaPlotLabel = "genus")

#----------------------------------------------------#
# Replot ANCOM-BC results
#----------------------------------------------------#

source("00-functions.R") # for ancomReplot()
#--------------WIS--------------#
# HMS
ancomReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            fontSize = 7,
            pointSize = 0.5)
ancomReplot(seqCenter = "HMS",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_WIS",
            fontSize = 7,
            pointSize = 0.5)
# BCM
ancomReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS")
ancomReplot(seqCenter = "BCM",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_WIS")
# WashU
ancomReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "WashU",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_WIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# Broad_WGS
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            fontSize = 7,
            pointSize = 0.5)
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_WIS",
            fontSize = 7,
            pointSize = 0.5)
# MDA
ancomReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "MDA",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_WIS",
            numticks = 2,
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# UNC
ancomReplot(seqCenter = "UNC",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            plotWidth = 12,
            plotHeight = 5,
            fontSize = 6,
            nRowFlag = TRUE,
            nRow = 3,
            pointSize = 0.5)
# CMS
ancomReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)

#--------------BinsWIS--------------#
# HMS
ancomReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            fontSize = 7,
            pointSize = 0.5)
ancomReplot(seqCenter = "HMS",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_BinsWIS",
            fontSize = 7,
            pointSize = 0.5)
# BCM
ancomReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS")
ancomReplot(seqCenter = "BCM",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_BinsWIS")
# WashU
ancomReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "WashU",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_BinsWIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# Broad_WGS
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            fontSize = 7,
            pointSize = 0.5)
ancomReplot(seqCenter = "Broad_WGS",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_BinsWIS",
            fontSize = 7,
            pointSize = 0.5)
# MDA
ancomReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
ancomReplot(seqCenter = "MDA",
            sampleType = "Blood Derived Normal",
            dataString = "K1C_HPRC_BinsWIS",
            numticks = 2,
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
# UNC
ancomReplot(seqCenter = "UNC",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            plotWidth = 12,
            plotHeight = 5,
            fontSize = 6,
            nRowFlag = TRUE,
            nRow = 3,
            pointSize = 0.5)
# CMS
ancomReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            plotWidth = 6,
            fontSize = 6,
            pointSize = 0.5)
#----------------------------------------------------#
# Replot alpha diversity results
#----------------------------------------------------#

source("00-functions.R") # for alphaReplot()
#--------------WIS--------------#
# HMS
alphaReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# BCM
alphaReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            alphaDivType = "Observed",
            plotWidth = 6,
            plotHeight = 3,
            fontSize = 8)
# Broad_WGS
alphaReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# WashU
alphaReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            alphaDivType = "Observed",
            plotWidth = 6,
            plotHeight = 3,
            fontSize = 8)
# MDA
alphaReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# CMS
alphaReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_WIS",
            alphaDivType = "Observed",
            plotWidth = 4,
            plotHeight = 3,
            fontSize = 8,
            ptOnlyFlag = TRUE)

#--------------BinsWIS--------------#
# HMS
alphaReplot(seqCenter = "HMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# BCM
alphaReplot(seqCenter = "BCM",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            alphaDivType = "Observed",
            plotWidth = 6,
            plotHeight = 3,
            fontSize = 8)
# Broad_WGS
alphaReplot(seqCenter = "Broad_WGS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# WashU
alphaReplot(seqCenter = "WashU",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            alphaDivType = "Observed",
            plotWidth = 6,
            plotHeight = 3,
            fontSize = 8)
# MDA
alphaReplot(seqCenter = "MDA",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            alphaDivType = "Observed",
            plotWidth = 8,
            plotHeight = 3,
            fontSize = 8)
# CMS
alphaReplot(seqCenter = "CMS",
            sampleType = "Primary Tumor",
            dataString = "K1C_HPRC_BinsWIS",
            alphaDivType = "Observed",
            plotWidth = 4,
            plotHeight = 3,
            fontSize = 8,
            ptOnlyFlag = TRUE)

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



