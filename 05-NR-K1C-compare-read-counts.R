# 05-NR-K1C-compare-read-counts.R
# Author: Greg Poore
# Date: Feb 23, 2024
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
require(ggpubr)
require(ggsci)

numCores <- detectCores()
registerDoMC(cores=numCores)

# Create figure folder
plotPrefix <- "K1C-read-cmp"
if(!( dir.exists( file.path( paste0("Figures/",plotPrefix) )))){
  dir.create(file.path( paste0("Figures/",plotPrefix) ))
}
if(!( dir.exists( file.path( paste0("Interim_data/",plotPrefix) )))){
  dir.create(file.path( paste0("Interim_data/",plotPrefix) ))
}

#----------------------------------------------------#
# Import data
#----------------------------------------------------#

load("Input_data/tcgaVbDataAndMetadataAndSNM.RData", verbose = TRUE)

load("Interim_data/K1C-HPRC/data_vsnm_tcga_full_wis_bins_features_subset_20Feb24.RData",verbose = TRUE)

countsVbFinalNonzeroQC_HPRC <- countsVbFinalNonzeroQC
countsVbFinalNonzeroQCWIS_HPRC <- countsVbFinalNonzeroQCWIS
metaDataFinalNonzeroQC_HPRC <- metaDataFinalNonzeroQC

rm(countsVbFinalNonzeroQC,
   voomDataGenusKrakenQCFilt,
   vsnmDataGenusKrakenQCFilt,
   countsVbFinalNonzeroQCWIS,
   voomDataGenusKrakenQCFiltWIS,
   vsnmDataGenusKrakenQCFiltWIS,
   countsVbFinalNonzeroQCBins,
   voomDataGenusKrakenQCFiltBins,
   vsnmDataGenusKrakenQCFiltBins,
   metaDataFinalNonzeroQC)

#----------------------------------------------------#
# Import data
#----------------------------------------------------#
load("Interim_data/K1C-T2T/data_raw_tcga_full_wis_bins_features_subset_23Feb24.RData",verbose = TRUE)

countsVbFinalNonzeroQC_T2T <- countsVbFinalNonzeroQC
countsVbFinalNonzeroQCWIS_T2T <- countsVbFinalNonzeroQCWIS
metaDataFinalNonzeroQC_T2T <- metaDataFinalNonzeroQC

rm(countsVbFinalNonzeroQC,
   countsVbFinalNonzeroQCWIS,
   countsVbFinalNonzeroQCBins,
   metaDataFinalNonzeroQC)

#----------------------------------------------------#
# Deriving single table of reads per sample
#----------------------------------------------------#
# Orig data
perSampleReads_hg19 <- data.frame(knightID = rownames(vbDataBarnDFReconciledQC),
                                  sample_type = metadataSamplesAllQC$sample_type,
                                  experimental_strategy = factor(metadataSamplesAllQC$experimental_strategy,
                                                  levels = c("WGS","RNA-Seq")),
                             investigation = gsub("^TCGA-","",metadataSamplesAllQC$investigation),
                             hg19 = rowSums(vbDataBarnDFReconciledQC))
# Format T2T data
countsVbFinalNonzeroQC_T2T_Form <- countsVbFinalNonzeroQC_T2T
rownames(countsVbFinalNonzeroQC_T2T_Form) <- metaDataFinalNonzeroQC_T2T$knightlabID
perSampleReads_T2T <- data.frame(knightID = rownames(countsVbFinalNonzeroQC_T2T_Form),
                                 T2T = rowSums(countsVbFinalNonzeroQC_T2T_Form))

# Format HPRC data
countsVbFinalNonzeroQC_HPRC_Form <- countsVbFinalNonzeroQC_HPRC
rownames(countsVbFinalNonzeroQC_HPRC_Form) <- metaDataFinalNonzeroQC_HPRC$knightlabID
perSampleReads_HPRC <- data.frame(knightID = rownames(countsVbFinalNonzeroQC_HPRC_Form),
                                 HPRC = rowSums(countsVbFinalNonzeroQC_HPRC_Form))

sum(rownames(countsVbFinalNonzeroQC_HPRC_Form) %in% perSampleReads_hg19$knightID) # 15180

# Left join
perSampleReads_All <- perSampleReads_hg19 %>%
  left_join(perSampleReads_T2T, by = "knightID") %>%
  left_join(perSampleReads_HPRC, by = "knightID") %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("knightID")

perSampleReads_All  %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal", "Solid Tissue Normal")) %>%
  mutate(sample_type = factor(case_when(
    sample_type == "Primary Tumor" ~ "PT",
    sample_type == "Blood Derived Normal" ~ "BDN",
    sample_type == "Solid Tissue Normal" ~ "STN"
  ), levels = c("PT","STN","BDN"))) %>%
  reshape2::melt(id.vars = c("experimental_strategy","sample_type","investigation")) %>%
  ggboxplot(x = "sample_type",
            y = "value",
            fill = "variable",
            outlier.size = 0.5,
            palette = c("#0072B5FF","#E18727FF","#BC3C29FF"),
            ylab = "Per-sample microbial counts (Kraken 1)",
            xlab = "",
            legend = "top",
            # add = "mean",
            facet.by = "experimental_strategy") +
  labs(fill = "Sequential host depletion:") +
  scale_y_log10(limits = c(1e1,10^9.5), breaks=c(1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9))
ggsave(filename = paste0("Figures/",plotPrefix,"/readDepthChange_hg19_vs_T2T_vs_Pan_23Oct24.jpeg"),
       width = 6, height = 6, units = "in", dpi = "retina")



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
