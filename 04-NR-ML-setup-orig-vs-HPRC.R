# 04-NR-ML-setup-orig-vs-HPRC.R
# Author: Greg Poore
# Date: Feb 20, 2024
# Purpose: Compare ML models built with 2020 Kraken data vs 2024 Kraken data

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
plotPrefix <- "Orig-vs-HPRC"
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
load("Input_data/snmDataSampleTypeWithExpStrategyFINAL.RData", verbose = TRUE)
metadataSamplesAllQC$cgc_platform <- metadataSamplesAllQC$platform

# Load HPRC data
load("Interim_data/K1C-HPRC/data_vsnm_tcga_full_wis_bins_features_subset_20Feb24.RData", verbose=TRUE)
# countsVbFinalNonzeroQC
# voomDataGenusKrakenQCFilt
# vsnmDataGenusKrakenQCFilt
# countsVbFinalNonzeroQCWIS
# voomDataGenusKrakenQCFiltWIS
# vsnmDataGenusKrakenQCFiltWIS
# countsVbFinalNonzeroQCBins
# voomDataGenusKrakenQCFiltBins
# vsnmDataGenusKrakenQCFiltBins
# metaDataFinalNonzeroQC

load("Interim_data/K1C-HPRC/data_per_center_tcga_full_wis_bins_features_subset_20Feb24.RData", verbose=TRUE)

#----------------------------------------------------------#
# Identify overlaps and subset
#----------------------------------------------------------#
# Make rownames of new data = knight lab IDs
all(rownames(countsVbFinalNonzeroQC) == rownames(metaDataFinalNonzeroQC)) # TRUE
countsVbFinalNonzeroQCForm <- countsVbFinalNonzeroQC
rownames(countsVbFinalNonzeroQCForm) <- metaDataFinalNonzeroQC$knightlabID

dim(countsVbFinalNonzeroQCForm) # 15180  1874
dim(vbDataBarnDFReconciledQC) # 17625  1993

sum(colnames(countsVbFinalNonzeroQCForm) %in% colnames(vbDataBarnDFReconciledQC)) # 1873 (ie, 120 diff from 1993)
sum(colnames(vbDataBarnDFReconciledQC) %in% colnames(countsVbFinalNonzeroQCForm)) # 1873 (ie, 120 diff from 1993)
#--> so only 1 feature is new in countsVbFinalNonzeroQCForm
colnames(countsVbFinalNonzeroQCForm)[which(!(colnames(countsVbFinalNonzeroQCForm) %in% colnames(vbDataBarnDFReconciledQC)))]
#--> new one is k__Viruses.o__Picornavirales.f__Picornaviridae.g__Teschovirus

## Save feature sets for later intersection and enrichment
totalFeatureSetRaw <- unique(c(colnames(vbDataBarnDFReconciledQC),
                               colnames(countsVbFinalNonzeroQCForm)))
length(totalFeatureSetRaw) # 1994
totalFeatureSetVSNM <- unique(c(colnames(snmDataSampleTypeWithExpStrategy),
                                colnames(vsnmDataGenusKrakenQCFilt)))
length(totalFeatureSetVSNM) # 1891

# save(totalFeatureSetRaw,
#      totalFeatureSetVSNM,
#      file = "Interim_data/totalFeatureSet_old_vs_new_21Feb24.RData")

## Subset orig meta, count, and VSNM data to same samples
metadataSamplesAllQC_MatchedSamp <- droplevels(metadataSamplesAllQC[rownames(countsVbFinalNonzeroQCForm),])
vbDataBarnDFReconciledQC_MatchedSamp <- vbDataBarnDFReconciledQC[rownames(countsVbFinalNonzeroQCForm),]
dim(vbDataBarnDFReconciledQC_MatchedSamp) # 15180  1993
snmDataSampleTypeWithExpStrategy_MatchedSamp <- snmDataSampleTypeWithExpStrategy[rownames(countsVbFinalNonzeroQCForm),]
dim(snmDataSampleTypeWithExpStrategy_MatchedSamp) # 15180 1795

# Subset new data (count and VSNM) by features
countsVbFinalNonzeroQC_MatchedFeat <- countsVbFinalNonzeroQC[,colnames(countsVbFinalNonzeroQC) %in%
                                                               colnames(snmDataSampleTypeWithExpStrategy_MatchedSamp)]
dim(countsVbFinalNonzeroQC_MatchedFeat) # 15180  1778

vsnmDataGenusKrakenQCFilt_MatchedFeat <- vsnmDataGenusKrakenQCFilt[,colnames(vsnmDataGenusKrakenQCFilt) %in%
                                                                     colnames(snmDataSampleTypeWithExpStrategy_MatchedSamp)]
dim(vsnmDataGenusKrakenQCFilt_MatchedFeat) # 15180  1778

## Outline comparisons
# 1) VSNM orig 17625x1795 vs VSNM New 15180x1874 --> ready
# 2) VSNM orig 17625x1795 vs VSNM New 15180x1778 --> ready
# 3) VSNM orig 15180x1795 vs VSNM New 15180x1874 --> ready
# 4) VSNM orig 15180x1795 vs VSNM New 15180x1778 --> ready
# 5) VSNM orig re-norm 15180x1993 vs VSNM New 15180x1874 --> below
# 6) VSNM orig re-norm 15180x1993 vs VSNM New re-norm 15180x1778 --> below
# 7) VSNM orig 15180x1795 vs VSNM New re-norm 15180x1778 --> below
# 8) Per center: Raw orig 17625x1993 vs Raw New 14914x1874 --> below
# 9) Per center: Raw orig 14914x1993 vs Raw New 14914x1874 --> below

#-------------------------------------------#
# (5) VSNM Re-norm orig data 15180x1795
#-------------------------------------------#

source("00-functions.R")
# Note: covDesignNorm had dimensions: 15180    18
batchCorrectedVSNMDataRenormOrig <- vsnmFunctionTCGA(qcData = vbDataBarnDFReconciledQC_MatchedSamp,
                                                     qcMetadata = metadataSamplesAllQC_MatchedSamp)
snmDataSampleTypeWithExpStrategy_RenormMatchedSamp <- batchCorrectedVSNMDataRenormOrig$snmData

#-------------------------------------------#
# (6-7) VSNM Re-norm new data 15180x1778
#-------------------------------------------#

source("00-functions.R")
# Note: covDesignNorm had dimensions: 17625   211
batchCorrectedVSNMDataRenormNew <- vsnmFunctionTCGA(qcData = countsVbFinalNonzeroQC_MatchedFeat,
                                                    qcMetadata = metaDataFinalNonzeroQC)
vsnmDataGenusKrakenQCFilt_RenormMatchedFeat <- batchCorrectedVSNMDataRenormNew$snmData

#-------------------------------------------#
# (8) Per-center data for Raw orig 17625x1993
#-------------------------------------------#

#-----------------Subset to Illumina HiSeq-----------------#
metadataSamplesAllQC_HiSeq <- metadataSamplesAllQC %>%
  filter(platform == "Illumina HiSeq") %>% droplevels()

vbDataBarnDFReconciledQC_HiSeq <- vbDataBarnDFReconciledQC[rownames(metadataSamplesAllQC_HiSeq),]
#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_HiSeq_HMS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_BCM <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_MDA <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_WashU <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_Broad_WGS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metadataSamplesAllQC_HiSeq_UNC <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metadataSamplesAllQC_HiSeq_CMS <- metadataSamplesAllQC_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

# WGS
metadataSamplesAllQC_HiSeq_WGS <- metadataSamplesAllQC_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA
metadataSamplesAllQC_HiSeq_RNA <- metadataSamplesAllQC_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset Full count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
vbDataBarnDFReconciledQC_HiSeq_HMS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_HMS),]
vbDataBarnDFReconciledQC_HiSeq_BCM <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_BCM),]
vbDataBarnDFReconciledQC_HiSeq_MDA <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_MDA),]
vbDataBarnDFReconciledQC_HiSeq_WashU <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WashU),]
vbDataBarnDFReconciledQC_HiSeq_Broad_WGS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_Broad_WGS),]
vbDataBarnDFReconciledQC_HiSeq_WGS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
vbDataBarnDFReconciledQC_HiSeq_UNC <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_UNC),]
vbDataBarnDFReconciledQC_HiSeq_CMS <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_CMS),]
vbDataBarnDFReconciledQC_HiSeq_RNA <- vbDataBarnDFReconciledQC_HiSeq[rownames(metadataSamplesAllQC_HiSeq_RNA),]

#-------------------------------------------#
# (9) Per-center data for Raw orig 14914x1993
#-------------------------------------------#

#-----------------Subset to Illumina HiSeq-----------------#
metadataSamplesAllQC_MatchedSamp_HiSeq <- metadataSamplesAllQC %>%
  rownames_to_column("knightID") %>%
  filter(knightID %in% metaDataFinalNonzeroQC_HiSeq$knightlabID) %>%
  column_to_rownames("knightID") %>% droplevels()

vbDataBarnDFReconciledQC_MatchedSamp_HiSeq <- vbDataBarnDFReconciledQC[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq),]
#--------------------Subset metadata by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metadataSamplesAllQC_MatchedSamp_HiSeq_HMS <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(data_submitting_center_label == "Harvard Medical School") %>% 
  droplevels()
metadataSamplesAllQC_MatchedSamp_HiSeq_BCM <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% 
  droplevels()
metadataSamplesAllQC_MatchedSamp_HiSeq_MDA <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% 
  droplevels()
metadataSamplesAllQC_MatchedSamp_HiSeq_WashU <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% 
  droplevels()
metadataSamplesAllQC_MatchedSamp_HiSeq_Broad_WGS <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% 
  droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, but RNA-Seq is only GBM, only a WGS is made)
metadataSamplesAllQC_MatchedSamp_HiSeq_UNC <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(data_submitting_center_label == "University of North Carolina") %>% 
  droplevels()
metadataSamplesAllQC_MatchedSamp_HiSeq_CMS <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% 
  droplevels()

# WGS
metadataSamplesAllQC_MatchedSamp_HiSeq_WGS <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA
metadataSamplesAllQC_MatchedSamp_HiSeq_RNA <- metadataSamplesAllQC_MatchedSamp_HiSeq %>% 
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset Full count data by seq center--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_HMS <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_HMS),]
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_BCM <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_BCM),]
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_MDA <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_MDA),]
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_WashU <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_WashU),]
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_Broad_WGS <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_Broad_WGS),]
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_WGS <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_UNC <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_UNC),]
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_CMS <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_CMS),]
vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_RNA <- vbDataBarnDFReconciledQC_MatchedSamp_HiSeq[rownames(metadataSamplesAllQC_MatchedSamp_HiSeq_RNA),]

## Outline comparisons
# 1) VSNM orig 17625x1795 vs VSNM New 15180x1874 --> ready
# 2) VSNM orig 17625x1795 vs VSNM New 15180x1778 --> ready
# 3) VSNM orig 15180x1795 vs VSNM New 15180x1874 --> ready
# 4) VSNM orig 15180x1795 vs VSNM New 15180x1778 --> ready
# 5) VSNM orig re-norm 15180x1993 vs VSNM New 15180x1874 --> below
# 6) VSNM orig re-norm 15180x1993 vs VSNM New re-norm 15180x1778 --> below
# 7) VSNM orig 15180x1795 vs VSNM New re-norm 15180x1778 --> below
# 8) Per center: Raw orig 17625x1993 vs Raw New 14914x1874 --> below
# 9) Per center: Raw orig 14914x1993 vs Raw New 14914x1874 --> below

save(
  #----------Old----------#
  
  # VSNM orig 17625x1795
  snmDataSampleTypeWithExpStrategy,
  
  # VSNM orig 15180x1795
  snmDataSampleTypeWithExpStrategy_MatchedSamp,
  
  # VSNM orig re-norm 15180x1993
  snmDataSampleTypeWithExpStrategy_RenormMatchedSamp,
  
  # Per center: Raw orig 17625x1993
  vbDataBarnDFReconciledQC_HiSeq_HMS,
  vbDataBarnDFReconciledQC_HiSeq_BCM,
  vbDataBarnDFReconciledQC_HiSeq_MDA,
  vbDataBarnDFReconciledQC_HiSeq_WashU,
  vbDataBarnDFReconciledQC_HiSeq_Broad_WGS,
  vbDataBarnDFReconciledQC_HiSeq_UNC,
  vbDataBarnDFReconciledQC_HiSeq_CMS,
  
  vbDataBarnDFReconciledQC_HiSeq_WGS,
  vbDataBarnDFReconciledQC_HiSeq_RNA,
  
  # Per center: Raw orig 14914x1993
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_HMS,
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_BCM,
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_MDA,
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_WashU,
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_Broad_WGS,
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_UNC,
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_CMS,
  
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_WGS,
  vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_RNA,
  
  # Subset metadata
  metadataSamplesAllQC_HiSeq_HMS,
  metadataSamplesAllQC_HiSeq_BCM,
  metadataSamplesAllQC_HiSeq_MDA,
  metadataSamplesAllQC_HiSeq_WashU,
  metadataSamplesAllQC_HiSeq_Broad_WGS,
  metadataSamplesAllQC_HiSeq_UNC,
  metadataSamplesAllQC_HiSeq_CMS,
  
  # Subset metadata
  metadataSamplesAllQC_MatchedSamp_HiSeq_HMS,
  metadataSamplesAllQC_MatchedSamp_HiSeq_BCM,
  metadataSamplesAllQC_MatchedSamp_HiSeq_MDA,
  metadataSamplesAllQC_MatchedSamp_HiSeq_WashU,
  metadataSamplesAllQC_MatchedSamp_HiSeq_Broad_WGS,
  metadataSamplesAllQC_MatchedSamp_HiSeq_UNC,
  metadataSamplesAllQC_MatchedSamp_HiSeq_CMS,
  
  # WGS and RNA subset metadata
  metadataSamplesAllQC_HiSeq_WGS,
  metadataSamplesAllQC_HiSeq_RNA,
  metadataSamplesAllQC_MatchedSamp_HiSeq_WGS,
  metadataSamplesAllQC_MatchedSamp_HiSeq_RNA,
  
  # Full metadata
  metadataSamplesAllQC,
  metadataSamplesAllQC_MatchedSamp,
  metadataSamplesAllQC_HiSeq,
  
  #----------New----------#
  
  # VSNM New 15180x1874
  vsnmDataGenusKrakenQCFilt,
  
  # VSNM New 15180x1778
  vsnmDataGenusKrakenQCFilt_MatchedFeat,
  
  # VSNM New re-norm 15180x1778
  vsnmDataGenusKrakenQCFilt_RenormMatchedFeat,
  
  # Raw New 14914x1874
  countsVbFinalNonzeroQC_HiSeq_HMS,
  countsVbFinalNonzeroQC_HiSeq_BCM,
  countsVbFinalNonzeroQC_HiSeq_MDA,
  countsVbFinalNonzeroQC_HiSeq_WashU,
  countsVbFinalNonzeroQC_HiSeq_Broad_WGS,
  countsVbFinalNonzeroQC_HiSeq_UNC,
  countsVbFinalNonzeroQC_HiSeq_CMS,
  
  countsVbFinalNonzeroQC_HiSeq_WGS,
  countsVbFinalNonzeroQC_HiSeq_RNA,
  
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
  metaDataFinalNonzeroQC,
  metaDataFinalNonzeroQC_HiSeq,
  file = paste0("Interim_data/",plotPrefix,"/data_tcga_all_comparisons_old_vs_new_20Feb24.RData"))

## Outline comparisons
# 1) VSNM orig 17625x1795 vs VSNM New 15180x1874 --> ready
# 2) VSNM orig 17625x1795 vs VSNM New 15180x1778 --> ready
# 3) VSNM orig 15180x1795 vs VSNM New 15180x1874 --> ready
# 4) VSNM orig 15180x1795 vs VSNM New 15180x1778 --> ready
# 5) VSNM orig re-norm 15180x1993 vs VSNM New 15180x1874 --> below
# 6) VSNM orig re-norm 15180x1993 vs VSNM New re-norm 15180x1778 --> below
# 7) VSNM orig 15180x1795 vs VSNM New re-norm 15180x1778 --> below
# 8) Per center: Raw orig 17625x1993 vs Raw New 14914x1874 --> below
# 9) Per center: Raw orig 14914x1993 vs Raw New 14914x1874 --> below

#-------------------------------------------#
#-------------------------------------------#
# Feature enrichment and similarity
#-------------------------------------------#
#-------------------------------------------#

load("Interim_data/totalFeatureSet_old_vs_new_21Feb24.RData", verbose = TRUE)
abbreviationsTCGA_Allcancer <- read.csv("Input_data/tcga_abbreviations.csv", 
                                        stringsAsFactors = FALSE, row.names = 1)

## Write function
runEnrichmentFxn3_BatchCorrected <- function(dataSetOldInput,
                                             dataSetNewInput,
                                             fileNameStringInput="OldvsNew1-VSNM-orig-17625x1795-vs-VSNM-New-15180x1874",
                                             totalFeaturesInput=totalFeatureSetVSNM){
  source("00-functions.R") # for enrichmentFxn2() function
  # PT
  print("Working on PT...")
  ovn_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                           totalFeatures = totalFeaturesInput,
                           seqCenter = "",
                           sampleType = "Primary Tumor",
                           dataSetRaw = dataSetOldInput, # OLD
                           dataSetVSNM = dataSetNewInput, # NEW
                           cancerAbbrevs = abbreviationsTCGA_Allcancer,
                           plotWidth = 8,
                           fileNameString = fileNameStringInput,
                           kendallOnlyFlag = FALSE,
                           showCM = FALSE)
  
  # BDN
  print("Working on BDN...")
  ovn_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                            totalFeatures = totalFeaturesInput,
                            seqCenter = "",
                            sampleType = "Blood Derived Normal",
                            dataSetRaw = dataSetOldInput, # OLD
                            dataSetVSNM = dataSetNewInput, # NEW
                            cancerAbbrevs = abbreviationsTCGA_Allcancer,
                            plotWidth = 6,
                            fileNameString = fileNameStringInput,
                            kendallOnlyFlag = FALSE,
                            showCM = FALSE)
  
  # PT vs STN
  print("Working on PT vs STN...")
  ovn_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                            totalFeatures = totalFeaturesInput,
                            seqCenter = "",
                            sampleType = "Primary Tumor vs Solid Tissue Normal",
                            dataSetRaw = dataSetOldInput, # OLD
                            dataSetVSNM = dataSetNewInput, # NEW
                            cancerAbbrevs = abbreviationsTCGA_Allcancer,
                            plotWidth = 4,
                            fileNameString = fileNameStringInput,
                            kendallOnlyFlag = FALSE,
                            showCM = FALSE)
  
  print("Overall results:")
  resPT <- ovn_PT$fisherKendallCombinedDf
  resBDN <- ovn_BDN$fisherKendallCombinedDf
  resSTN <- ovn_STN$fisherKendallCombinedDf
  
  print(sprintf("PT Fisher: %d / %d significant overlap",
                sum(resPT$p.adj.signif.fisher != "ns"), nrow(resPT)))
  print(sprintf("PT Kendall: %d / %d significant overlap",
                sum(resPT$p.adj.signif.kendall != "ns"), nrow(resPT)))
  
  print(sprintf("BDN Fisher: %d / %d significant overlap",
                sum(resBDN$p.adj.signif.fisher != "ns"), nrow(resBDN)))
  print(sprintf("BDN Kendall: %d / %d significant overlap",
                sum(resBDN$p.adj.signif.kendall != "ns"), nrow(resBDN)))
  
  print(sprintf("STN Fisher: %d / %d significant overlap",
                sum(resSTN$p.adj.signif.fisher != "ns"), nrow(resSTN)))
  print(sprintf("STN Kendall: %d / %d significant overlap",
                sum(resSTN$p.adj.signif.kendall != "ns"), nrow(resSTN)))
  
  res <- list(ovn_PT=ovn_PT,
              ovn_BDN=ovn_BDN,
              ovn_STN=ovn_STN)
  return(res)
}

#-------------------------------------------#
# 1) VSNM orig 17625x1795 vs VSNM New 15180x1874
# snmDataSampleTypeWithExpStrategy vs vsnmDataGenusKrakenQCFilt
#-------------------------------------------#
source("00-functions.R")
cmp1Res <- runEnrichmentFxn3_BatchCorrected(dataSetOldInput="snmDataSampleTypeWithExpStrategy",
                                 dataSetNewInput="vsnmDataGenusKrakenQCFilt",
                                 fileNameStringInput="OldvsNew1-VSNM-orig-17625x1795-vs-VSNM-New-15180x1874",
                                 totalFeaturesInput=totalFeatureSetVSNM)
# [1] "Overall results:"
# [1] "PT Fisher: 28 / 32 significant overlap"
# [1] "PT Kendall: 29 / 32 significant overlap"
# [1] "BDN Fisher: 14 / 20 significant overlap"
# [1] "BDN Kendall: 14 / 20 significant overlap"
# [1] "STN Fisher: 6 / 14 significant overlap"
# [1] "STN Kendall: 6 / 14 significant overlap"

#-------------------------------------------#
# 2) VSNM orig 17625x1795 vs VSNM New 15180x1778
# snmDataSampleTypeWithExpStrategy vs vsnmDataGenusKrakenQCFilt_MatchedFeat
#-------------------------------------------#

source("00-functions.R")
cmp2Res <- runEnrichmentFxn3_BatchCorrected(dataSetOldInput="snmDataSampleTypeWithExpStrategy",
                                 dataSetNewInput="vsnmDataGenusKrakenQCFilt_MatchedFeat",
                                 fileNameStringInput="OldvsNew2-VSNM-orig-17625x1795-vs-VSNM-New-15180x1778",
                                 totalFeaturesInput=totalFeatureSetVSNM)
# [1] "Overall results:"
# [1] "PT Fisher: 27 / 32 significant overlap"
# [1] "PT Kendall: 29 / 32 significant overlap"
# [1] "BDN Fisher: 15 / 20 significant overlap"
# [1] "BDN Kendall: 15 / 20 significant overlap"
# [1] "STN Fisher: 6 / 14 significant overlap"
# [1] "STN Kendall: 7 / 14 significant overlap"

#-------------------------------------------#
# 3) VSNM orig 15180x1795 vs VSNM New 15180x1874
# snmDataSampleTypeWithExpStrategy_MatchedSamp vs vsnmDataGenusKrakenQCFilt
#-------------------------------------------#

source("00-functions.R")
cmp3Res <- runEnrichmentFxn3_BatchCorrected(dataSetOldInput="snmDataSampleTypeWithExpStrategy_MatchedSamp",
                                 dataSetNewInput="vsnmDataGenusKrakenQCFilt",
                                 fileNameStringInput="OldvsNew3-VSNM-orig-15180x1795-vs-VSNM-New-15180x1874",
                                 totalFeaturesInput=totalFeatureSetVSNM)
# [1] "Overall results:"
# [1] "PT Fisher: 29 / 32 significant overlap"
# [1] "PT Kendall: 29 / 32 significant overlap"
# [1] "BDN Fisher: 16 / 20 significant overlap"
# [1] "BDN Kendall: 16 / 20 significant overlap"
# [1] "STN Fisher: 7 / 14 significant overlap"
# [1] "STN Kendall: 7 / 14 significant overlap"

#-------------------------------------------#
# 4) VSNM orig 15180x1795 vs VSNM New 15180x1778
# snmDataSampleTypeWithExpStrategy_MatchedSamp vs vsnmDataGenusKrakenQCFilt_MatchedFeat
#-------------------------------------------#

source("00-functions.R")
cmp4Res <- runEnrichmentFxn3_BatchCorrected(dataSetOldInput="snmDataSampleTypeWithExpStrategy_MatchedSamp",
                                            dataSetNewInput="vsnmDataGenusKrakenQCFilt_MatchedFeat",
                                            fileNameStringInput="OldvsNew4-VSNM-orig-15180x1795-vs-VSNM-New-15180x1778",
                                            totalFeaturesInput=totalFeatureSetVSNM)
# [1] "Overall results:"
# [1] "PT Fisher: 28 / 32 significant overlap"
# [1] "PT Kendall: 28 / 32 significant overlap"
# [1] "BDN Fisher: 16 / 20 significant overlap"
# [1] "BDN Kendall: 17 / 20 significant overlap"
# [1] "STN Fisher: 7 / 14 significant overlap"
# [1] "STN Kendall: 7 / 14 significant overlap"

#-------------------------------------------#
# 5) VSNM orig re-norm 15180x1993 vs VSNM New 15180x1874
# snmDataSampleTypeWithExpStrategy_RenormMatchedSamp vs vsnmDataGenusKrakenQCFilt
#-------------------------------------------#

source("00-functions.R")
cmp5Res <- runEnrichmentFxn3_BatchCorrected(dataSetOldInput="snmDataSampleTypeWithExpStrategy_RenormMatchedSamp",
                                            dataSetNewInput="vsnmDataGenusKrakenQCFilt",
                                            fileNameStringInput="OldvsNew5-VSNM-orig-15180x1993-vs-VSNM-New-15180x1874",
                                            totalFeaturesInput=totalFeatureSetRaw)
# [1] "Overall results:"
# [1] "PT Fisher: 23 / 32 significant overlap"
# [1] "PT Kendall: 24 / 32 significant overlap"
# [1] "BDN Fisher: 16 / 20 significant overlap"
# [1] "BDN Kendall: 16 / 20 significant overlap"
# [1] "STN Fisher: 7 / 14 significant overlap"
# [1] "STN Kendall: 9 / 14 significant overlap"

#-------------------------------------------#
# 6) VSNM orig re-norm 15180x1993 vs VSNM New re-norm 15180x1778
# snmDataSampleTypeWithExpStrategy_RenormMatchedSamp vs vsnmDataGenusKrakenQCFilt_RenormMatchedFeat
#-------------------------------------------#

source("00-functions.R")
cmp6Res <- runEnrichmentFxn3_BatchCorrected(dataSetOldInput="snmDataSampleTypeWithExpStrategy_RenormMatchedSamp",
                                            dataSetNewInput="vsnmDataGenusKrakenQCFilt_RenormMatchedFeat",
                                            fileNameStringInput="OldvsNew6-VSNM-orig-15180x1993-vs-VSNM-New-15180x1778",
                                            totalFeaturesInput=totalFeatureSetRaw)
# [1] "Overall results:"
# [1] "PT Fisher: 23 / 32 significant overlap"
# [1] "PT Kendall: 23 / 32 significant overlap"
# [1] "BDN Fisher: 16 / 20 significant overlap"
# [1] "BDN Kendall: 16 / 20 significant overlap"
# [1] "STN Fisher: 0 / 14 significant overlap"
# [1] "STN Kendall: 3 / 14 significant overlap"

#-------------------------------------------#
# 7) VSNM orig 15180x1795 vs VSNM New re-norm 15180x1778
# snmDataSampleTypeWithExpStrategy_MatchedSamp vs vsnmDataGenusKrakenQCFilt_RenormMatchedFeat
#-------------------------------------------#

source("00-functions.R")
cmp7Res <- runEnrichmentFxn3_BatchCorrected(dataSetOldInput="snmDataSampleTypeWithExpStrategy_MatchedSamp",
                                            dataSetNewInput="vsnmDataGenusKrakenQCFilt_RenormMatchedFeat",
                                            fileNameStringInput="OldvsNew7-VSNM-orig-15180x1795-vs-VSNM-New-15180x1778",
                                            totalFeaturesInput=totalFeatureSetVSNM)
# [1] "Overall results:"
# [1] "PT Fisher: 29 / 32 significant overlap"
# [1] "PT Kendall: 30 / 32 significant overlap"
# [1] "BDN Fisher: 16 / 20 significant overlap"
# [1] "BDN Kendall: 16 / 20 significant overlap"
# [1] "STN Fisher: 6 / 14 significant overlap"
# [1] "STN Kendall: 7 / 14 significant overlap"

#-------------------------------------------#
# 8) Per center: Raw orig 17625x1993 vs Raw New 14914x1874
# vbDataBarnDFReconciledQC_HiSeq_HMS,
# vbDataBarnDFReconciledQC_HiSeq_BCM,
# vbDataBarnDFReconciledQC_HiSeq_MDA,
# vbDataBarnDFReconciledQC_HiSeq_WashU,
# vbDataBarnDFReconciledQC_HiSeq_Broad_WGS,
# vbDataBarnDFReconciledQC_HiSeq_UNC,
# vbDataBarnDFReconciledQC_HiSeq_CMS

# countsVbFinalNonzeroQC_HiSeq_HMS,
# countsVbFinalNonzeroQC_HiSeq_BCM,
# countsVbFinalNonzeroQC_HiSeq_MDA,
# countsVbFinalNonzeroQC_HiSeq_WashU,
# countsVbFinalNonzeroQC_HiSeq_Broad_WGS,
# countsVbFinalNonzeroQC_HiSeq_UNC,
# countsVbFinalNonzeroQC_HiSeq_CMS
#-------------------------------------------#
source("00-functions.R") # for enrichmentFxn2() function

#-----------------HMS-----------------#
# PT
ovn8_HMS_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                         totalFeatures = totalFeatureSetRaw,
                         seqCenter = "HMS",
                         sampleType = "Primary Tumor",
                         dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                         dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                         plotWidth = 5,
                         fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# BDN
ovn8_HMS_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                          totalFeatures = totalFeatureSetRaw,
                          seqCenter = "HMS",
                          sampleType = "Blood Derived Normal",
                          dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                          dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                          plotWidth = 5,
                          fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn8_HMS_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                          totalFeatures = totalFeatureSetRaw,
                          seqCenter = "HMS",
                          sampleType = "Primary Tumor vs Solid Tissue Normal",
                          dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                          dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                          plotWidth = 2.5,
                          fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

#-----------------BCM-----------------#

# PT
ovn8_BCM_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "BCM",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 4,
                             fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# BDN
ovn8_BCM_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "BCM",
                              sampleType = "Blood Derived Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 4,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn8_BCM_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "BCM",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 3,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

#-----------------MDA-----------------#
# PT
ovn8_MDA_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "MDA",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 4,
                             fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# BDN
ovn8_MDA_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "MDA",
                              sampleType = "Blood Derived Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 4,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

#-----------------WashU-----------------#
# PT
ovn8_WashU_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "WashU",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 3,
                             fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# BDN
ovn8_WashU_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "WashU",
                              sampleType = "Blood Derived Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 3,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

#-----------------Broad_WGS-----------------#
# PT
ovn8_Broad_WGS_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "Broad_WGS",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 5,
                             fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# BDN
ovn8_Broad_WGS_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "Broad_WGS",
                              sampleType = "Blood Derived Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 4,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn8_Broad_WGS_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "Broad_WGS",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 2.5,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

#-----------------UNC-----------------#
# PT
ovn8_UNC_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "UNC",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 10,
                             fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn8_UNC_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "UNC",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 6,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

#-----------------CMS-----------------#
# PT
ovn8_CMS_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "CMS",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 3,
                             fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn8_CMS_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "CMS",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 2.5,
                              fileNameString = "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874")

#----------------Combine outputs----------------#
kf_Comb <- rbind(# HMS
  ovn8_HMS_BDN$fisherKendallCombinedDf,
  ovn8_HMS_PT$fisherKendallCombinedDf,
  ovn8_HMS_STN$fisherKendallCombinedDf,
  # BCM
  ovn8_BCM_BDN$fisherKendallCombinedDf,
  ovn8_BCM_PT$fisherKendallCombinedDf,
  ovn8_BCM_STN$fisherKendallCombinedDf,
  # MDA
  ovn8_MDA_BDN$fisherKendallCombinedDf,
  ovn8_MDA_PT$fisherKendallCombinedDf,
  # WashU
  ovn8_WashU_BDN$fisherKendallCombinedDf,
  ovn8_WashU_PT$fisherKendallCombinedDf,
  # Broad_WGS
  ovn8_Broad_WGS_BDN$fisherKendallCombinedDf,
  ovn8_Broad_WGS_PT$fisherKendallCombinedDf,
  ovn8_Broad_WGS_STN$fisherKendallCombinedDf,
  # UNC
  ovn8_UNC_PT$fisherKendallCombinedDf,
  ovn8_UNC_STN$fisherKendallCombinedDf,
  # CMS
  ovn8_CMS_PT$fisherKendallCombinedDf,
  ovn8_CMS_STN$fisherKendallCombinedDf)

kf_CombP <- kf_Comb %>%
  group_by(abbrev, ST) %>%
  summarise(p.fisher.comb = survcomp::combine.test(p.fisher),
            p.kendall.comb = survcomp::combine.test(p.kendall),
            tau.comb = median(tau),
            tau.se = ifelse(is.na(sd(tau)/n()),0,sd(tau)/n()),
            OR.comb = median(OR),
            OR.se = ifelse(is.na(sd(OR)/n()),0,sd(OR)/n()),
            count = n()) %>%
  rstatix::adjust_pvalue("p.fisher.comb", "p.fisher.comb.adj") %>%
  rstatix::adjust_pvalue("p.kendall.comb", "p.kendall.comb.adj") %>%
  rstatix::add_significance(p.col = "p.fisher.comb.adj") %>%
  rstatix::add_significance(p.col = "p.kendall.comb.adj")

fileNameString <- "OldvsNew8-Raw-orig-17625x1993-vs-Raw-New-14914x1874"
seqCenter <- "All"
plotWidthInput <- c(8,6,4)
sampleTypeInput <- c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal")
for(ii in 1:length(sampleTypeInput)){
  sampleType <- sampleTypeInput[ii]
  plotWidth <- plotWidthInput[ii]
  kf_CombP_Filt <- kf_CombP %>%
    filter(ST == sampleType)
  
  p.adj.signif.fisher.kendall <- c(kf_CombP_Filt$p.fisher.comb.adj.signif,
                                   kf_CombP_Filt$p.kendall.comb.adj.signif)
  p.adj.signif.kendall.fisher <- c(kf_CombP_Filt$p.kendall.comb.adj.signif,
                                   kf_CombP_Filt$p.fisher.comb.adj.signif)
  
  kf_CombP_Filt %>%
    mutate(log.p.fisher.comb.adj = -log10(p.fisher.comb.adj),
           log.p.kendall.comb.adj = -log10(p.kendall.comb.adj)) %>%
    select(abbrev, log.p.fisher.comb.adj, log.p.kendall.comb.adj,
           p.fisher.comb.adj.signif, p.kendall.comb.adj.signif, count) %>%
    reshape2::melt(id.vars = c("abbrev","p.fisher.comb.adj.signif","p.kendall.comb.adj.signif","count")) %>%
    mutate(variable = factor(case_when(
      variable == "log.p.fisher.comb.adj" ~ "Fisher",
      variable == "log.p.kendall.comb.adj" ~ "Kendall",
    ), levels = c("Kendall","Fisher"))) %>%
    mutate(p.adj.signif.combined = `if`(head(as.character(variable),1)=="Fisher",
                                        p.adj.signif.fisher.kendall,
                                        p.adj.signif.kendall.fisher)) %>%
    ggplot(aes(x = reorder(abbrev,-value,median), y = value, fill=variable)) +
    geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
    geom_text(aes(label = p.adj.signif.combined, y = value),
              position = position_dodge(0.9), 
              hjust = -0.2,
              vjust = 0.8,
              size = 3,
              angle = 90) +
    geom_text(aes(label = count, y = value),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
    ylim(c(0,1.1*max( c(-log10(kf_CombP_Filt$p.fisher.comb.adj),
                        -log10(kf_CombP_Filt$p.kendall.comb.adj)) ))) +
    scale_fill_nejm(name = "Test type") +
    theme_pubr() +
    rotate_x_text(90) +
    theme(axis.text.x = element_text(size=10)) +
    labs(x = "TCGA Cancer Type", 
         y = "-Log(p-adjust)",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> pvalPlot
  
  print(pvalPlot)
  fileNamePval <- paste0("Figures/",fileNameString,"/pvalue_combined_Fisher_Kendall_barplot_tcga_",seqCenter,"_",
                         gsub('([[:punct:]])|\\s+','',sampleType),
                         ".jpeg")
  ggsave(filename = fileNamePval,
         plot = pvalPlot,
         dpi = "retina", units = "in", height = 3, width = plotWidth)
  
  kf_CombP_Filt %>%
    ggplot(aes(x = reorder(abbrev, -tau.comb,median), y = tau.comb)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=tau.comb-tau.se, ymax=tau.comb+tau.se), width=.2) +
    theme_pubr() +
    geom_text(aes(label = count, y = tau.comb),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    rotate_x_text(90) +
    labs(x = "TCGA Cancer Type", 
         y = "Kendall tau",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    geom_text(aes(label = p.kendall.comb.adj.signif, y = tau.comb), vjust = -0.4) +
    ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb) ))) +
    # ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb+kf_CombP_Filt$tau.se) ))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotKendallTau
  
  print(barPlotKendallTau)
  fileNameTau <- paste0("Figures/",fileNameString,"/kendall_combined_tau_barplot_tcga_",seqCenter,"_",
                        gsub('([[:punct:]])|\\s+','',sampleType),
                        ".jpeg")
  ggsave(filename = fileNameTau,
         plot = barPlotKendallTau,
         dpi = "retina", units = "in", height = 3.5, width = plotWidth)
  
  kf_CombP_Filt %>%
    ggplot(aes(x = reorder(abbrev, -OR.comb,median), y = OR.comb)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=OR.comb-OR.se, ymax=OR.comb+OR.se), width=.2) +
    theme_pubr() +
    geom_text(aes(label = count, y = OR.comb),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    rotate_x_text(90) +
    labs(x = "TCGA Cancer Type", 
         y = "Odds ratio feature enrichment",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    geom_text(aes(label = p.fisher.comb.adj.signif, y = OR.comb), vjust = -0.4) +
    ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb) ))) +
    # ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb+kf_CombP_Filt$OR.se) ))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotOR
  
  print(barPlotOR)
  fileNameOR <- paste0("Figures/",fileNameString,"/odds_ratio_combined_barplot_tcga_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
  ggsave(filename = fileNameOR,
         plot = barPlotOR,
         dpi = "retina", units = "in", height = 3.5, width = plotWidth)
}

#-------------------------------------------#
# 9) Per center: Raw orig 14914x1993 vs Raw New 14914x1874
# vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_HMS,
# vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_BCM,
# vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_MDA,
# vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_WashU,
# vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_Broad_WGS,
# vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_UNC,
# vbDataBarnDFReconciledQC_MatchedSamp_HiSeq_CMS,

# countsVbFinalNonzeroQC_HiSeq_HMS,
# countsVbFinalNonzeroQC_HiSeq_BCM,
# countsVbFinalNonzeroQC_HiSeq_MDA,
# countsVbFinalNonzeroQC_HiSeq_WashU,
# countsVbFinalNonzeroQC_HiSeq_Broad_WGS,
# countsVbFinalNonzeroQC_HiSeq_UNC,
# countsVbFinalNonzeroQC_HiSeq_CMS
#-------------------------------------------#
source("00-functions.R") # for enrichmentFxn2() function

#-----------------HMS-----------------#
# PT
ovn9_HMS_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "HMS",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 5,
                             fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# BDN
ovn9_HMS_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "HMS",
                              sampleType = "Blood Derived Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 5,
                              fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn9_HMS_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "HMS",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 2.5,
                              fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

#-----------------BCM-----------------#

# PT
ovn9_BCM_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "BCM",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 4,
                             fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# BDN
ovn9_BCM_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "BCM",
                              sampleType = "Blood Derived Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 4,
                              fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn9_BCM_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "BCM",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 3,
                              fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

#-----------------MDA-----------------#
# PT
ovn9_MDA_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "MDA",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 4,
                             fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# BDN
ovn9_MDA_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "MDA",
                              sampleType = "Blood Derived Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 4,
                              fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

#-----------------WashU-----------------#
# PT
ovn9_WashU_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                               totalFeatures = totalFeatureSetRaw,
                               seqCenter = "WashU",
                               sampleType = "Primary Tumor",
                               dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                               dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                               plotWidth = 3,
                               fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# BDN
ovn9_WashU_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                                totalFeatures = totalFeatureSetRaw,
                                seqCenter = "WashU",
                                sampleType = "Blood Derived Normal",
                                dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                                dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                                plotWidth = 3,
                                fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

#-----------------Broad_WGS-----------------#
# PT
ovn9_Broad_WGS_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                                   totalFeatures = totalFeatureSetRaw,
                                   seqCenter = "Broad_WGS",
                                   sampleType = "Primary Tumor",
                                   dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                                   dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                                   plotWidth = 5,
                                   fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# BDN
ovn9_Broad_WGS_BDN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                                    totalFeatures = totalFeatureSetRaw,
                                    seqCenter = "Broad_WGS",
                                    sampleType = "Blood Derived Normal",
                                    dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                                    dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                                    plotWidth = 4,
                                    fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn9_Broad_WGS_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                                    totalFeatures = totalFeatureSetRaw,
                                    seqCenter = "Broad_WGS",
                                    sampleType = "Primary Tumor vs Solid Tissue Normal",
                                    dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                                    dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                                    plotWidth = 2.5,
                                    fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

#-----------------UNC-----------------#
# PT
ovn9_UNC_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "UNC",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 10,
                             fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn9_UNC_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "UNC",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 6,
                              fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

#-----------------CMS-----------------#
# PT
ovn9_CMS_PT <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                             totalFeatures = totalFeatureSetRaw,
                             seqCenter = "CMS",
                             sampleType = "Primary Tumor",
                             dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                             dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                             plotWidth = 3,
                             fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

# PT vs STN
ovn9_CMS_STN <- enrichmentFxn3(myPath = "Supporting_scripts/S05-ML-old-vs-new/features__",
                              totalFeatures = totalFeatureSetRaw,
                              seqCenter = "CMS",
                              sampleType = "Primary Tumor vs Solid Tissue Normal",
                              dataSetRaw = "vbDataBarnDFReconciledQC_MatchedSamp_HiSeq", # OLD
                              dataSetVSNM = "countsVbFinalNonzeroQC_HiSeq", # NEW
                              plotWidth = 2.5,
                              fileNameString = "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874")

#----------------Combine outputs----------------#
kf_Comb <- rbind(# HMS
  ovn9_HMS_BDN$fisherKendallCombinedDf,
  ovn9_HMS_PT$fisherKendallCombinedDf,
  ovn9_HMS_STN$fisherKendallCombinedDf,
  # BCM
  ovn9_BCM_BDN$fisherKendallCombinedDf,
  ovn9_BCM_PT$fisherKendallCombinedDf,
  ovn9_BCM_STN$fisherKendallCombinedDf,
  # MDA
  ovn9_MDA_BDN$fisherKendallCombinedDf,
  ovn9_MDA_PT$fisherKendallCombinedDf,
  # WashU
  ovn9_WashU_BDN$fisherKendallCombinedDf,
  ovn9_WashU_PT$fisherKendallCombinedDf,
  # Broad_WGS
  ovn9_Broad_WGS_BDN$fisherKendallCombinedDf,
  ovn9_Broad_WGS_PT$fisherKendallCombinedDf,
  ovn9_Broad_WGS_STN$fisherKendallCombinedDf,
  # UNC
  ovn9_UNC_PT$fisherKendallCombinedDf,
  ovn9_UNC_STN$fisherKendallCombinedDf,
  # CMS
  ovn9_CMS_PT$fisherKendallCombinedDf,
  ovn9_CMS_STN$fisherKendallCombinedDf)

kf_CombP <- kf_Comb %>%
  group_by(abbrev, ST) %>%
  summarise(p.fisher.comb = survcomp::combine.test(p.fisher),
            p.kendall.comb = survcomp::combine.test(p.kendall),
            tau.comb = median(tau),
            tau.se = ifelse(is.na(sd(tau)/n()),0,sd(tau)/n()),
            OR.comb = median(OR),
            OR.se = ifelse(is.na(sd(OR)/n()),0,sd(OR)/n()),
            count = n()) %>%
  rstatix::adjust_pvalue("p.fisher.comb", "p.fisher.comb.adj") %>%
  rstatix::adjust_pvalue("p.kendall.comb", "p.kendall.comb.adj") %>%
  rstatix::add_significance(p.col = "p.fisher.comb.adj") %>%
  rstatix::add_significance(p.col = "p.kendall.comb.adj")

fileNameString <- "OldvsNew9-Raw-orig-14914x1993-vs-Raw-New-14914x1874"
seqCenter <- "All"
plotWidthInput <- c(10,8,5)
sampleTypeInput <- c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal")
for(ii in 1:length(sampleTypeInput)){
  sampleType <- sampleTypeInput[ii]
  plotWidth <- plotWidthInput[ii]
  kf_CombP_Filt <- kf_CombP %>%
    filter(ST == sampleType)
  
  p.adj.signif.fisher.kendall <- c(kf_CombP_Filt$p.fisher.comb.adj.signif,
                                   kf_CombP_Filt$p.kendall.comb.adj.signif)
  p.adj.signif.kendall.fisher <- c(kf_CombP_Filt$p.kendall.comb.adj.signif,
                                   kf_CombP_Filt$p.fisher.comb.adj.signif)
  
  kf_CombP_Filt %>%
    mutate(log.p.fisher.comb.adj = -log10(p.fisher.comb.adj),
           log.p.kendall.comb.adj = -log10(p.kendall.comb.adj)) %>%
    select(abbrev, log.p.fisher.comb.adj, log.p.kendall.comb.adj,
           p.fisher.comb.adj.signif, p.kendall.comb.adj.signif, count) %>%
    reshape2::melt(id.vars = c("abbrev","p.fisher.comb.adj.signif","p.kendall.comb.adj.signif","count")) %>%
    mutate(variable = factor(case_when(
      variable == "log.p.fisher.comb.adj" ~ "Fisher",
      variable == "log.p.kendall.comb.adj" ~ "Kendall",
    ), levels = c("Kendall","Fisher"))) %>%
    mutate(p.adj.signif.combined = `if`(head(as.character(variable),1)=="Fisher",
                                        p.adj.signif.fisher.kendall,
                                        p.adj.signif.kendall.fisher)) %>%
    ggplot(aes(x = reorder(abbrev,-value,median), y = value, fill=variable)) +
    geom_bar(stat="identity", color="black", position=position_dodge(0.9)) +
    geom_text(aes(label = p.adj.signif.combined, y = value),
              position = position_dodge(0.9), 
              hjust = -0.2,
              vjust = 0.8,
              size = 3,
              angle = 90) +
    geom_text(aes(label = count, y = value),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
    ylim(c(0,1.1*max( c(-log10(kf_CombP_Filt$p.fisher.comb.adj),
                        -log10(kf_CombP_Filt$p.kendall.comb.adj)) ))) +
    scale_fill_nejm(name = "Test type") +
    theme_pubr() +
    rotate_x_text(90) +
    theme(axis.text.x = element_text(size=10)) +
    labs(x = "TCGA Cancer Type", 
         y = "-Log(p-adjust)",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> pvalPlot
  
  print(pvalPlot)
  fileNamePval <- paste0("Figures/",fileNameString,"/pvalue_combined_Fisher_Kendall_barplot_tcga_",seqCenter,"_",
                         gsub('([[:punct:]])|\\s+','',sampleType),
                         ".jpeg")
  ggsave(filename = fileNamePval,
         plot = pvalPlot,
         dpi = "retina", units = "in", height = 3, width = plotWidth)
  
  kf_CombP_Filt %>%
    ggplot(aes(x = reorder(abbrev, -tau.comb,median), y = tau.comb)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=tau.comb-tau.se, ymax=tau.comb+tau.se), width=.2) +
    theme_pubr() +
    geom_text(aes(label = count, y = tau.comb),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    rotate_x_text(90) +
    labs(x = "TCGA Cancer Type", 
         y = "Kendall tau",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    geom_text(aes(label = p.kendall.comb.adj.signif, y = tau.comb), vjust = -0.4) +
    ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb) ))) +
    # ylim(c(0,1.1*max( (kf_CombP_Filt$tau.comb+kf_CombP_Filt$tau.se) ))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotKendallTau
  
  print(barPlotKendallTau)
  fileNameTau <- paste0("Figures/",fileNameString,"/kendall_combined_tau_barplot_tcga_",seqCenter,"_",
                        gsub('([[:punct:]])|\\s+','',sampleType),
                        ".jpeg")
  ggsave(filename = fileNameTau,
         plot = barPlotKendallTau,
         dpi = "retina", units = "in", height = 3.5, width = plotWidth)
  
  kf_CombP_Filt %>%
    ggplot(aes(x = reorder(abbrev, -OR.comb,median), y = OR.comb)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=OR.comb-OR.se, ymax=OR.comb+OR.se), width=.2) +
    theme_pubr() +
    geom_text(aes(label = count, y = OR.comb),
              position = position_dodge(0.9), 
              hjust = 0.5,
              vjust = 1.5,
              size = 2,
              angle = 0,
              color = "white") +
    rotate_x_text(90) +
    labs(x = "TCGA Cancer Type", 
         y = "Odds ratio feature enrichment",
         title = paste(seqCenter,sampleType,sep = " | ")) +
    geom_text(aes(label = p.fisher.comb.adj.signif, y = OR.comb), vjust = -0.4) +
    ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb) ))) +
    # ylim(c(0,1.1*max( (kf_CombP_Filt$OR.comb+kf_CombP_Filt$OR.se) ))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right") -> barPlotOR
  
  print(barPlotOR)
  fileNameOR <- paste0("Figures/",fileNameString,"/odds_ratio_combined_barplot_tcga_",seqCenter,"_",
                       gsub('([[:punct:]])|\\s+','',sampleType),
                       ".jpeg")
  ggsave(filename = fileNameOR,
         plot = barPlotOR,
         dpi = "retina", units = "in", height = 3.5, width = plotWidth)
}

#-------------------------------------------#
#-------------------------------------------#
# ML plots
#-------------------------------------------#
#-------------------------------------------#
# Going with #1 and #8 above

#-------------------------------------------#
# 1) VSNM orig 17625x1795 vs VSNM New 15180x1874
# snmDataSampleTypeWithExpStrategy vs vsnmDataGenusKrakenQCFilt
#-------------------------------------------#

## Load data
mlPerf_OldvsNew_Cmp1 <- read.csv("Supporting_scripts/S05-ML-old-vs-new/perfML_Old_vs_New_ALL_20Feb24.csv", 
                                               stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=aucroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5) %>%
  filter(datasetName %in% c("snmDataSampleTypeWithExpStrategy","vsnmDataGenusKrakenQCFilt"))

# AUROC
mlPerf_OldvsNew_Cmp1 %>%
  select(abbrev,datasetName,sampleType,AUROC) %>%
  tidyr::spread(key = datasetName, value = AUROC) %>%
  filter(complete.cases(.)) %>%
  ggscatter(x = "snmDataSampleTypeWithExpStrategy",
            y = "vsnmDataGenusKrakenQCFilt",
            xlab = "Original AUROC (VSNM data)",
            ylab = "Updated AUROC (VSNM data)",
            conf.int = TRUE,
            conf.int.level = 0.95,
            color = "sampleType",
            palette = "nejm"
            ) +
  stat_smooth(method = "lm", se = T, fullrange = T, color = "black") +
  stat_cor(method = "spearman",cor.coef.name = "rho")
ggsave(filename = paste0("Figures/",plotPrefix,"/oldVsNew_Corr_Cmp1_AUROC_21Feb24_spearman.jpeg"),
       dpi = "retina", units = "in", width = 3, height = 4) # width used to be 4.5

# AUPR
mlPerf_OldvsNew_Cmp1 %>%
  select(abbrev,datasetName,sampleType,AUPR) %>%
  mutate(sampleType = gsub("Primary Tumor vs Solid Tissue Normal","Tumor vs Normal",sampleType)) %>%
  tidyr::spread(key = datasetName, value = AUPR) %>%
  filter(complete.cases(.)) %>%
  ggscatter(x = "snmDataSampleTypeWithExpStrategy",
            y = "vsnmDataGenusKrakenQCFilt",
            xlab = "Original AUPR (VSNM data)",
            ylab = "Updated AUPR (VSNM data)",
            conf.int = TRUE,
            ylim = c(0,1),
            conf.int.level = 0.95,
            color = "sampleType",
            palette = "nejm"
  ) +
  stat_smooth(method = "lm", se = T, fullrange = T, color = "black") +
  # theme(legend.position = "right") +
  stat_cor(
    method = "spearman",cor.coef.name = "rho")
ggsave(filename = paste0("Figures/",plotPrefix,"/oldVsNew_Corr_Cmp1_AUPR_21Feb24_spearman.jpeg"),
       dpi = "retina", units = "in", width = 3, height = 4) # width used to be 4.5

#-------------------------------------------#
# 8) Per center: Raw orig 17625x1993 vs Raw New 14914x1874
# vbDataBarnDFReconciledQC_HiSeq_HMS,
# vbDataBarnDFReconciledQC_HiSeq_BCM,
# vbDataBarnDFReconciledQC_HiSeq_MDA,
# vbDataBarnDFReconciledQC_HiSeq_WashU,
# vbDataBarnDFReconciledQC_HiSeq_Broad_WGS,
# vbDataBarnDFReconciledQC_HiSeq_UNC,
# vbDataBarnDFReconciledQC_HiSeq_CMS

# countsVbFinalNonzeroQC_HiSeq_HMS,
# countsVbFinalNonzeroQC_HiSeq_BCM,
# countsVbFinalNonzeroQC_HiSeq_MDA,
# countsVbFinalNonzeroQC_HiSeq_WashU,
# countsVbFinalNonzeroQC_HiSeq_Broad_WGS,
# countsVbFinalNonzeroQC_HiSeq_UNC,
# countsVbFinalNonzeroQC_HiSeq_CMS
#-------------------------------------------#
source("Supporting_scripts/S00-SummarySE.R")

## Load data
mlPerf_OldvsNew_Cmp8 <- read.csv("Supporting_scripts/S05-ML-old-vs-new/perfML_Old_vs_New_ALL_20Feb24.csv", 
                                 stringsAsFactors = FALSE) %>%
  select(-X) %>%
  distinct() %>%
  rename(AUROC=aucroc, AUPR=aupr) %>%
  mutate(abbrev = abbreviationsTCGA_Allcancer[diseaseType,"abbrev"]) %>%
  mutate(nullAUPR = ifelse(minorityClassName == "SolidTissueNormal",
                           yes=majorityClassSize/(minorityClassSize+majorityClassSize),
                           no=minorityClassSize/(minorityClassSize+majorityClassSize))) %>%
  mutate(nullAUROC = 0.5) %>%
  filter(grepl("vbDataBarnDFReconciledQC_HiSeq_|countsVbFinalNonzeroQC_HiSeq_",datasetName)) %>%
  mutate(seqCenter = gsub("^vbDataBarnDFReconciledQC_HiSeq_|^countsVbFinalNonzeroQC_HiSeq_","",datasetName))

mlPerf_OldvsNew_Cmp8 %>%
  mutate(datasetName = gsub("_HMS|_BCM|_MDA|_WashU|_Broad_WGS|_UNC|_CMS","",datasetName)) %>%
  mutate(datasetName = gsub("countsVbFinalNonzeroQC_HiSeq","Updated",datasetName)) %>%
  mutate(datasetName = gsub("vbDataBarnDFReconciledQC_HiSeq","Original",datasetName)) %>%
  mutate(sampleType = gsub("Primary Tumor vs Solid Tissue Normal","Tumor vs Normal",sampleType)) %>%
  select(AUROC, AUPR, abbrev,diseaseType,
         sampleType,datasetName, seqCenter) %>%
  ggboxplot(x = "seqCenter",
            y = "AUROC",
            fill = "datasetName",
            palette = c("#0072B5FF","#BC3C29FF")) +
  facet_grid( ~ sampleType, scales = "free") +
  rotate_x_text(30) +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1), limits = c(0.4,1))
ggsave(filename = paste0("Figures/",plotPrefix,"/oldVsNew_Corr_Cmp8_AUROC_21Feb24.jpeg"),
       dpi = "retina", units = "in", width = 8, height = 5)

mlPerf_OldvsNew_Cmp8 %>%
  mutate(datasetName = gsub("_HMS|_BCM|_MDA|_WashU|_Broad_WGS|_UNC|_CMS","",datasetName)) %>%
  mutate(datasetName = gsub("countsVbFinalNonzeroQC_HiSeq","Updated",datasetName)) %>%
  mutate(datasetName = gsub("vbDataBarnDFReconciledQC_HiSeq","Original",datasetName)) %>%
  mutate(sampleType = gsub("Primary Tumor vs Solid Tissue Normal","Tumor vs Normal",sampleType)) %>%
  select(AUROC, AUPR, abbrev,diseaseType,
         sampleType,datasetName, seqCenter) %>%
  ggboxplot(x = "seqCenter",
            y = "AUPR",
            fill = "datasetName",
            palette = c("#0072B5FF","#BC3C29FF")) +
  facet_grid( ~ sampleType, scales = "free") +
  rotate_x_text(30) +
  theme(text=element_text(size=15)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1))
ggsave(filename = paste0("Figures/",plotPrefix,"/oldVsNew_Corr_Cmp8_AUPR_21Feb24.jpeg"),
       dpi = "retina", units = "in", width = 8, height = 5)


