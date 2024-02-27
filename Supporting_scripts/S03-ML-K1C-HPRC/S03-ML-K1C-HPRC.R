#-----------------------------------------------------------------------------
# S03-ML-K1C-HPRC.R
# Copyright (c) 2024--, Greg Poore
# Purpose: Run machine learning on all available TCGA cancer types by seq center
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
require(splitstackshape)
require(reshape2)
require(tidyr)
require(caret) # for model building
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning
require(cvAUC)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("../../Interim_data/K1C-HPRC/data_per_center_tcga_full_wis_bins_features_subset_20Feb24.RData")
load("../../Interim_data/K1C-HPRC/data_vsnm_tcga_full_wis_bins_features_subset_20Feb24.RData")
# load("../../Interim_data/K1C-HPRC/metadata_stage_noMut_subsets_20Feb24.RData")

# GBM HYPERPARAMETER SEARCH GRID BELOW (DEFAULT PER CARET PACKAGE)
kfoldGBMGrid <- data.frame(n.trees=150, interaction.depth=3,
                           shrinkage=0.1,
                           n.minobsinnode=1)

# Model setup -- MODIFY AS NEEDED; ALL THESE COMPARISONS WILL BE TESTED
sampleTypeList <- c("Blood Derived Normal",
                    "Primary Tumor vs Solid Tissue Normal",
                    "Primary Tumor")
# MODIFY AS NEEDED
datasetList <- list(
  # # BinsWIS data per-center
  # countsVbFinalNonzeroQCBinsWIS_HiSeq_HMS,
  # countsVbFinalNonzeroQCBinsWIS_HiSeq_BCM,
  # countsVbFinalNonzeroQCBinsWIS_HiSeq_MDA,
  # countsVbFinalNonzeroQCBinsWIS_HiSeq_WashU,
  # countsVbFinalNonzeroQCBinsWIS_HiSeq_Broad_WGS,
  # countsVbFinalNonzeroQCBinsWIS_HiSeq_UNC,
  # countsVbFinalNonzeroQCBinsWIS_HiSeq_CMS,
  # 
  # # Bins data per-center
  # countsVbFinalNonzeroQCBins_HiSeq_HMS,
  # countsVbFinalNonzeroQCBins_HiSeq_BCM,
  # countsVbFinalNonzeroQCBins_HiSeq_MDA,
  # countsVbFinalNonzeroQCBins_HiSeq_WashU,
  # countsVbFinalNonzeroQCBins_HiSeq_Broad_WGS,
  # countsVbFinalNonzeroQCBins_HiSeq_UNC,
  # countsVbFinalNonzeroQCBins_HiSeq_CMS,
  # 
  # # WIS data per-center
  # countsVbFinalNonzeroQCWIS_HiSeq_HMS,
  # countsVbFinalNonzeroQCWIS_HiSeq_BCM,
  # countsVbFinalNonzeroQCWIS_HiSeq_MDA,
  # countsVbFinalNonzeroQCWIS_HiSeq_WashU,
  # countsVbFinalNonzeroQCWIS_HiSeq_Broad_WGS,
  # countsVbFinalNonzeroQCWIS_HiSeq_UNC,
  # countsVbFinalNonzeroQCWIS_HiSeq_CMS,
  # 
  # # Full data per-center
  # countsVbFinalNonzeroQC_HiSeq_HMS,
  # countsVbFinalNonzeroQC_HiSeq_BCM,
  # countsVbFinalNonzeroQC_HiSeq_MDA,
  # countsVbFinalNonzeroQC_HiSeq_WashU,
  # countsVbFinalNonzeroQC_HiSeq_Broad_WGS,
  # countsVbFinalNonzeroQC_HiSeq_UNC,
  # countsVbFinalNonzeroQC_HiSeq_CMS,
  
  # VSNM data
  vsnmDataGenusKrakenQCFiltWIS,
  vsnmDataGenusKrakenQCFilt
  )

datasetListNames <- c(
  # # BinsWIS data per-center
  # "countsVbFinalNonzeroQCBinsWIS_HiSeq_HMS",
  # "countsVbFinalNonzeroQCBinsWIS_HiSeq_BCM",
  # "countsVbFinalNonzeroQCBinsWIS_HiSeq_MDA",
  # "countsVbFinalNonzeroQCBinsWIS_HiSeq_WashU",
  # "countsVbFinalNonzeroQCBinsWIS_HiSeq_Broad_WGS",
  # "countsVbFinalNonzeroQCBinsWIS_HiSeq_UNC",
  # "countsVbFinalNonzeroQCBinsWIS_HiSeq_CMS",
  # 
  # # Bins data per-center
  # "countsVbFinalNonzeroQCBins_HiSeq_HMS",
  # "countsVbFinalNonzeroQCBins_HiSeq_BCM",
  # "countsVbFinalNonzeroQCBins_HiSeq_MDA",
  # "countsVbFinalNonzeroQCBins_HiSeq_WashU",
  # "countsVbFinalNonzeroQCBins_HiSeq_Broad_WGS",
  # "countsVbFinalNonzeroQCBins_HiSeq_UNC",
  # "countsVbFinalNonzeroQCBins_HiSeq_CMS",
  # 
  # # WIS data per-center
  # "countsVbFinalNonzeroQCWIS_HiSeq_HMS",
  # "countsVbFinalNonzeroQCWIS_HiSeq_BCM",
  # "countsVbFinalNonzeroQCWIS_HiSeq_MDA",
  # "countsVbFinalNonzeroQCWIS_HiSeq_WashU",
  # "countsVbFinalNonzeroQCWIS_HiSeq_Broad_WGS",
  # "countsVbFinalNonzeroQCWIS_HiSeq_UNC",
  # "countsVbFinalNonzeroQCWIS_HiSeq_CMS",
  # 
  # # Full data per-center
  # "countsVbFinalNonzeroQC_HiSeq_HMS",
  # "countsVbFinalNonzeroQC_HiSeq_BCM",
  # "countsVbFinalNonzeroQC_HiSeq_MDA",
  # "countsVbFinalNonzeroQC_HiSeq_WashU",
  # "countsVbFinalNonzeroQC_HiSeq_Broad_WGS",
  # "countsVbFinalNonzeroQC_HiSeq_UNC",
  # "countsVbFinalNonzeroQC_HiSeq_CMS",
  
  # VSNM data
  "vsnmDataGenusKrakenQCFiltWIS",
  "vsnmDataGenusKrakenQCFilt"
  )

metadataList <- list(
  # # BinsWIS data per-center
  # metaDataFinalNonzeroQC_HiSeq_HMS,
  # metaDataFinalNonzeroQC_HiSeq_BCM,
  # metaDataFinalNonzeroQC_HiSeq_MDA,
  # metaDataFinalNonzeroQC_HiSeq_WashU,
  # metaDataFinalNonzeroQC_HiSeq_Broad_WGS,
  # metaDataFinalNonzeroQC_HiSeq_UNC,
  # metaDataFinalNonzeroQC_HiSeq_CMS,
  # 
  # # Bins data per-center
  # metaDataFinalNonzeroQC_HiSeq_HMS,
  # metaDataFinalNonzeroQC_HiSeq_BCM,
  # metaDataFinalNonzeroQC_HiSeq_MDA,
  # metaDataFinalNonzeroQC_HiSeq_WashU,
  # metaDataFinalNonzeroQC_HiSeq_Broad_WGS,
  # metaDataFinalNonzeroQC_HiSeq_UNC,
  # metaDataFinalNonzeroQC_HiSeq_CMS,
  # 
  # # WIS data per-center
  # metaDataFinalNonzeroQC_HiSeq_HMS,
  # metaDataFinalNonzeroQC_HiSeq_BCM,
  # metaDataFinalNonzeroQC_HiSeq_MDA,
  # metaDataFinalNonzeroQC_HiSeq_WashU,
  # metaDataFinalNonzeroQC_HiSeq_Broad_WGS,
  # metaDataFinalNonzeroQC_HiSeq_UNC,
  # metaDataFinalNonzeroQC_HiSeq_CMS,
  # 
  # # Full data per-center
  # metaDataFinalNonzeroQC_HiSeq_HMS,
  # metaDataFinalNonzeroQC_HiSeq_BCM,
  # metaDataFinalNonzeroQC_HiSeq_MDA,
  # metaDataFinalNonzeroQC_HiSeq_WashU,
  # metaDataFinalNonzeroQC_HiSeq_Broad_WGS,
  # metaDataFinalNonzeroQC_HiSeq_UNC,
  # metaDataFinalNonzeroQC_HiSeq_CMS,
  
  # VSNM data
  metaDataFinalNonzeroQC,
  metaDataFinalNonzeroQC
  )

metadataListNames <- list(
  # # BinsWIS data per-center
  # "metaDataFinalNonzeroQC_HiSeq_HMS",
  # "metaDataFinalNonzeroQC_HiSeq_BCM",
  # "metaDataFinalNonzeroQC_HiSeq_MDA",
  # "metaDataFinalNonzeroQC_HiSeq_WashU",
  # "metaDataFinalNonzeroQC_HiSeq_Broad_WGS",
  # "metaDataFinalNonzeroQC_HiSeq_UNC",
  # "metaDataFinalNonzeroQC_HiSeq_CMS",
  # 
  # # Bins data per-center
  # "metaDataFinalNonzeroQC_HiSeq_HMS",
  # "metaDataFinalNonzeroQC_HiSeq_BCM",
  # "metaDataFinalNonzeroQC_HiSeq_MDA",
  # "metaDataFinalNonzeroQC_HiSeq_WashU",
  # "metaDataFinalNonzeroQC_HiSeq_Broad_WGS",
  # "metaDataFinalNonzeroQC_HiSeq_UNC",
  # "metaDataFinalNonzeroQC_HiSeq_CMS",
  # 
  # # WIS data per-center
  # "metaDataFinalNonzeroQC_HiSeq_HMS",
  # "metaDataFinalNonzeroQC_HiSeq_BCM",
  # "metaDataFinalNonzeroQC_HiSeq_MDA",
  # "metaDataFinalNonzeroQC_HiSeq_WashU",
  # "metaDataFinalNonzeroQC_HiSeq_Broad_WGS",
  # "metaDataFinalNonzeroQC_HiSeq_UNC",
  # "metaDataFinalNonzeroQC_HiSeq_CMS",
  # 
  # # Full data per-center
  # "metaDataFinalNonzeroQC_HiSeq_HMS",
  # "metaDataFinalNonzeroQC_HiSeq_BCM",
  # "metaDataFinalNonzeroQC_HiSeq_MDA",
  # "metaDataFinalNonzeroQC_HiSeq_WashU",
  # "metaDataFinalNonzeroQC_HiSeq_Broad_WGS",
  # "metaDataFinalNonzeroQC_HiSeq_UNC",
  # "metaDataFinalNonzeroQC_HiSeq_CMS",
  
  # VSNM data
  "metaDataFinalNonzeroQC",
  "metaDataFinalNonzeroQC"
  )

baseNamePerDatasetResultsFile <- "perfML_Human_HPRC_Full_WIS_Bins_20Feb24"
baseNameAllResultsFile <- "perfML_Human_HPRC_Full_WIS_Bins_ALL_20Feb24"

# MATCHED METADATA DATA FRAME FOR ALL COUNT DATA:
caretTuneGrid <- kfoldGBMGrid
numKFold <- 4
numResampleIter <- 1
prroc_roc <- list()
prroc_pr <- list()
perf <- list()
rep_perf <- list()
perfTmp <- list()
rep_perfTmp <- list()
perfTmp2 <- list()
rep_perfTmp2 <- list()

for(jj in seq_along(datasetList)){
  
  dataTmp <- datasetList[[jj]]
  datasetName <- datasetListNames[[jj]]
  
  metaTmpQC <- metadataList[[jj]]
  metaTmpQC$disease_type <- factor(metaTmpQC$disease_type)
  metadataName <- metadataListNames[[jj]]
  
  for(kk in seq_along(sampleTypeList)){
    st <- sampleTypeList[kk]
    print(st)
    
    if(st == "Primary Tumor vs Solid Tissue Normal"){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% c("Primary Tumor",
                                                                "Solid Tissue Normal"),])
    } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% st,])
    } else if(st == "Stage I vs IV"){
      metaTmp <- metaDataFinalNonzeroQCPath # metaTmpPath
      metaTmp2 <- droplevels(metaTmp[(metaTmp$sample_type %in% "Primary Tumor") &
                                       (metaTmp$pathologic_stage_label_binned %in% c("Stage1","Stage4")),])
      metaTmp2$disease_type <- factor(metaTmp2$disease_type)
    }
    
    print(seq_along(levels(metaTmp2$disease_type)))
    
    for(ii in seq_along(levels(metaTmp2$disease_type))){
      
      dt <- levels(metaTmp2$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- factor(gsub('([[:punct:]])|\\s+','',metaTmp3$sample_type))
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metaTmp3 <- metaTmp2
        metaTmp3$predY <- factor(ifelse(metaTmp2$disease_type == dt, 
                                        yes = dt, 
                                        no = "OtherCancerType"),
                                 levels = c(dt, "OtherCancerType"))
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      } else if(st == "Stage I vs IV"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- metaTmp3$pathologic_stage_label_binned
        positiveClass <- "Stage4"
        negativeClass <- "Stage1"
      }
      
      print(table(metaTmp3$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metaTmp3$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metaTmp3$predY) < 20)){next}
      
      minorityClassSize <- min(table((metaTmp3$predY)))
      majorityClassSize <- max(table((metaTmp3$predY)))
      
      minorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == min(table(metaTmp3$predY)))])
      majorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == max(table(metaTmp3$predY)))])
      
      mlDataY <- metaTmp3
      mlDataX <- dataTmp[rownames(mlDataY),]
      
      scrambleFlag <- ifelse(grepl("_scrambled",metadataName), yes = TRUE, no = FALSE)
      shuffleFlag <- ifelse(grepl("_shuffled",metadataName), yes = TRUE, no = FALSE)
      modelPreProcFlag <- ifelse(grepl("_VSNM",datasetName), yes = TRUE, no = FALSE)
      
      if(scrambleFlag){
        print("Scramble flag enabled! Scrambling metadata!")
        set.seed(42)
        mlDataY[,"predY"] <- sample(mlDataY[,"predY"])
      } else if(shuffleFlag){
        print("Shuffle flag enabled! Shuffling abundance data!")
        set.seed(42)
        mlDataX_shuffled <- mlDataX[sample(nrow(mlDataX)),]
        rownames(mlDataX_shuffled) <- rownames(mlDataX)
        mlDataX <- mlDataX_shuffled
      }
      
      set.seed(42)
      index <- createDataPartition(mlDataY$predY, p = 0.7, list = FALSE)
      trainX <- mlDataX[index,]
      trainY <- mlDataY[index,]$predY
      testX <- mlDataX[-index,]
      testY <- mlDataY[-index,]$predY
      refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
      refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
      # trainX <- mlDataX
      # trainY <- mlDataY[,"predY"]
      # refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           verboseIter = TRUE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      if(grepl("vsnm",datasetName,ignore.case = TRUE)){
        mlModel <- train(x = trainX,
                         y = refactoredTrainY,
                         method = "gbm",
                         preProcess = c("scale","center"),
                         trControl = ctrl,
                         verbose = TRUE,
                         metric = "ROC",
                         tuneGrid = caretTuneGrid)
      } else{
        mlModel <- train(x = trainX,
                         y = refactoredTrainY,
                         method = "gbm",
                         preProcess = c("nzv"),
                         trControl = ctrl,
                         verbose = TRUE,
                         metric = "ROC",
                         tuneGrid = caretTuneGrid)
      }
      
      
      predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
      fg <- predProbs[refactoredTestY == positiveClass]
      bg <- predProbs[refactoredTestY == negativeClass]
      
      prroc_roc[[ii]] <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      prroc_pr[[ii]] <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
      
      rocCurveData <- cbind(as.data.frame(prroc_roc[[ii]]$curve), disease_type = dt, sample_type = st)
      prCurveData <- cbind(as.data.frame(prroc_pr[[ii]]$curve), disease_type = dt, sample_type = st)
      
      perf[[ii]] <- data.frame(diseaseType = dt,
                               sampleType = st,
                               datasetName = datasetName,
                               metadataName = metadataName,
                               minorityClassSize = minorityClassSize,
                               majorityClassSize = majorityClassSize,
                               minorityClassName = minorityClassName,
                               majorityClassName = majorityClassName,
                               aucroc = prroc_roc[[ii]]$auc,
                               aupr = prroc_pr[[ii]]$auc.integral)
      
      print(perf[[ii]])
      
      confusionMatrix <- confusionMatrix(predict(mlModel, newdata = testX, type = "raw"), 
                                         refactoredTestY, 
                                         positive = positiveClass)
      
      print(confusionMatrix)
      
      # resPred <- mlModel$pred
      # 
      # ## Calculate performance on concatenated fold predictions
      # predProbs <- resPred
      # multiClass <- resPred
      # multiClass$pred <- relevel(multiClass$pred, positiveClass)
      # multiClass$obs <- relevel(multiClass$obs, positiveClass)
      # fg <- predProbs[resPred$obs == positiveClass,positiveClass]
      # bg <- predProbs[resPred$obs == negativeClass,positiveClass]
      # 
      # confusionMatrix <- confusionMatrix(multiClass$pred, multiClass$obs)
      # prroc_roc[[ii]] <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      # prroc_pr[[ii]] <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
      # rocCurveData <- cbind(as.data.frame(prroc_roc[[ii]]$curve), disease_type = dt, sample_type = st)
      # prCurveData <- cbind(as.data.frame(prroc_pr[[ii]]$curve), disease_type = dt, sample_type = st)
      # 
      # ## Estimate AUROC CIs using cvAUC on concatenated fold predictions
      # require(cvAUC)
      # out95 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.95)
      # out99 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.99)
      # # resCvAUC <- data.frame(estimate = c(out95$cvAUC, out99$cvAUC), se = c(out95$se, out99$se), lowerCI = c(out95$ci[1], out99$ci[1]), upperCI = c(out95$ci[1], out99$ci[2]), levelCI = c(0.95,0.99))
      # 
      # ## Split folds and calculate perf on each fold
      # resPredSplit <- split(resPred, resPred$Resample)
      # repX_perf <- list()
      # for(zz in seq_along(resPredSplit)){
      #   resPredSingleRep <- resPredSplit[[zz]]
      #   predProbs <- resPredSingleRep
      #   multiClass <- resPredSingleRep
      #   multiClass$pred <- relevel(multiClass$pred, positiveClass)
      #   multiClass$obs <- relevel(multiClass$obs, positiveClass)
      #   fg <- predProbs[resPredSingleRep$obs == positiveClass,positiveClass]
      #   bg <- predProbs[resPredSingleRep$obs == negativeClass,positiveClass]
      #   
      #   rep_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      #   rep_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
      #   
      #   repX_perf[[zz]] <- data.frame(auroc=rep_roc$auc,
      #                                 aupr=rep_pr$auc.integral,
      #                                 rep=paste0("Fold",zz), 
      #                                 diseaseType = dt,
      #                                 sampleType = st,
      #                                 datasetName = datasetName,
      #                                 metadataName = metadataName,
      #                                 minorityClassSize = minorityClassSize,
      #                                 majorityClassSize = majorityClassSize,
      #                                 minorityClassName = minorityClassName,
      #                                 majorityClassName = majorityClassName)
      # }
      # 
      # # SUMMARIZE MODEL PERFORMANCES
      # rep_perf[[ii]] <- do.call(rbind, repX_perf)
      # perf[[ii]] <- data.frame(auroc = prroc_roc[[ii]]$auc,
      #                          aupr = prroc_pr[[ii]]$auc.integral,
      #                          aucEstimate = out95$cvAUC, # either out95 or out99 work (same result)
      #                          aucSE95 = out95$se,
      #                          lowerCI95 = out95$ci[1],
      #                          upperCI95 = out95$ci[2],
      #                          aucSE99 = out99$se,
      #                          lowerCI99 = out99$ci[1],
      #                          upperCI99 = out99$ci[2],
      #                          diseaseType = dt,
      #                          sampleType = st,
      #                          datasetName = datasetName,
      #                          metadataName = metadataName,
      #                          minorityClassSize = minorityClassSize,
      #                          majorityClassSize = majorityClassSize,
      #                          minorityClassName = minorityClassName,
      #                          majorityClassName = majorityClassName)
      # 
      # print(perf[[ii]])
      # print(confusionMatrix)
      
      #--------------------------------------#
      # Save performance into relevant files #
      #--------------------------------------#
      
      filepathPerfPlots <- paste0("./perfPlots__",datasetName)
      filepathPerfPlotsDataPR <- paste0("./dataPR__",datasetName)
      filepathPerfPlotsDataROC <- paste0("./dataROC__",datasetName)
      filepathFeatures <- paste0("./features__",datasetName)
      filepathPerfStats <- paste0("./stats__",datasetName)
      filepathConfusionMatrix <- paste0("./confusionMatrix__",datasetName)
      
      if(!( dir.exists( file.path(filepathPerfPlots)))){
        dir.create(file.path(filepathPerfPlots))
      }
      if(!( dir.exists( file.path(filepathPerfPlotsDataPR)))){
        dir.create(file.path(filepathPerfPlotsDataPR))
      }
      if(!( dir.exists( file.path(filepathPerfPlotsDataROC)))){
        dir.create(file.path(filepathPerfPlotsDataROC))
      }
      if(!( dir.exists( file.path(filepathFeatures)))){
        dir.create(file.path(filepathFeatures))
      }
      if(!( dir.exists( file.path(filepathPerfStats)))){
        dir.create(file.path(filepathPerfStats))
      }
      if(!( dir.exists( file.path(filepathConfusionMatrix)))){
        dir.create(file.path(filepathConfusionMatrix))
      }
      
      filenameROC <- paste0(filepathPerfPlots,"/",
                            dt,
                            " -- ",
                            st,
                            " -- ROC.png")
      
      filenamePR <- paste0(filepathPerfPlots,"/",
                           dt,
                           " -- ",
                           st,
                           " -- PR.png")
      
      filenameROCData <- paste0(filepathPerfPlotsDataROC,"/",
                                dt,
                                " -- ",
                                st,
                                " -- ROC.csv")
      
      filenamePRData <- paste0(filepathPerfPlotsDataPR,"/",
                               dt,
                               " -- ",
                               st,
                               " -- PR.csv")
      
      filenameCM <- paste0(filepathConfusionMatrix,"/",
                           dt,
                           " -- ",
                           st,
                           " -- CM.txt")
      
      filenameFeatures <- paste0(filepathFeatures,"/",
                                 dt,
                                 " -- ",
                                 st,
                                 " -- Features.csv")
      
      filenameCSV <- paste0(filepathPerfStats,"/",
                            dt,
                            " -- ",
                            st,
                            " -- Perf.csv")
      
      write.csv(perf[[ii]], file = filenameCSV)
      
      png(filename=filenameROC, width = 6, height = 4, units = 'in', res = 300)
      plot(prroc_roc[[ii]])
      dev.off()
      
      png(filename=filenamePR, width = 6, height = 4, units = 'in', res = 300)
      plot(prroc_pr[[ii]])
      dev.off()
      
      write.table(prCurveData, sep=",", file = filenamePRData, col.names = FALSE)
      write.table(rocCurveData, sep=",", file = filenameROCData, col.names = FALSE)
      
      varImpBestModelDF <- as.data.frame(varImp( mlModel$finalModel, scale = FALSE ))
      varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
      varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
      colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
      varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
      write.csv(varImpBestModelDF2OrderedNonzero, file = filenameFeatures, row.names = FALSE)
      
      sink(file = filenameCM)
      print(dt)
      print(confusionMatrix)
      sink()
      
      rm(mlModel)
      
    }
    # SUMMARIZE RESULTS
    perfTmp[[kk]] <- do.call(rbind, perf)
  }
  # SUMMARIZE RESULTS
  perfTmp2[[jj]] <- do.call(rbind, perfTmp)
  
  write.csv(perfTmp2[[jj]], file = paste0(baseNamePerDatasetResultsFile,"_",datasetName,".csv"))
}

# SUMMARIZE RESULTS
perf1VsAll <- do.call(rbind, perfTmp2)

write.csv(perf1VsAll, file = paste0(baseNameAllResultsFile,".csv"))

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------