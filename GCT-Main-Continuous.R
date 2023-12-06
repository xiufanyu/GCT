################################################################################
##               Generalized Causal Tree (GCT) for Uplift Modeling            ## 
##                           (Continuous treatments)                          ##
################################################################################
rm(list=ls())
library(causalTree)
library(stringi)
library(data.tree)
library(stringr)
library(readr)
library(tidyr)
library(simcausal)
library(grf)
source('GCT-Utils-ConvertDataStructure.R')
source('GCT-Utils-OptimalCohorts.R')
source('GCT-Utils-Simulation-Continuous.R')

################################################################################
## Generate simulated data 
seed = 2
repId = 1
w1 = rnorm(1)
w2 = rnorm(1)
uncertaintyW <<- 1
suppressMessages(suppressWarnings(
  {
    trainingData <- generate_continuousT_training(seed = seed+repId, size = 500, 
                                uncertaintyW=uncertaintyW, w1, w2)
    trainingDataProcessed <- process_nonbinary_treatment(trainingData) 
    testingData <- generate_continuousT_training(seed = seed+repId+1, size = 500,
                               uncertaintyW=uncertaintyW, w1, w2)
    testingDataFeatures <- testingData[, c('ID', 'H1', 'H2', 'H3', 'H4')]
  }))
Avec  <- trainingData[, c('A')]
Abins <- rep(0, nrow(trainingData))
Abins[which(Avec <= 0.25 & Avec > 0)] <- 0.125
Abins[which(Avec <= 0.5 & Avec > 0.25)] <- 0.375
Abins[which(Avec <= 0.75 & Avec > 0.5)] <- 0.625
Abins[which(Avec <= 1 & Avec > 0.75)] <- 0.875
trainingDataBins <- data.frame(trainingData[,c('ID', 'H1', 'H2', 'H3', 'H4')], Abins, trainingData[,c('Y1')])
colnames(trainingDataBins) <- c('ID', 'H1', 'H2', 'H3', 'H4', 'A','Y1')
trainingDataProcessed <- process_nonbinary_treatment(trainingData)
trainingDataProcessedBins <- process_nonbinary_treatment(trainingDataBins)
trainingDataExport <- trainingDataBins[, c('H1', 'H2', 'A', 'Y1')]
colnames(trainingDataExport) <- c('H1', 'H2', 'A', 'Y1')
testingDataExport <- testingData[, c('H1', 'H2')]
################################################################################

################################################################################
## Task 1: Maximizing the average outcome with optimal treatment allocations
## over the testing data
  
## Method: GCT
tree_rpart_fit <- causalTree(Y1~H1+H2+H3+H4+Z,
                               data = trainingDataProcessed, 
                               treatment = trainingDataProcessed$treatment,
                               split.Rule = "CT", cv.option = "CT", 
                               split.Honest = T, cv.Honest = T, minsize = 20)
causalTreeDatatree <- convert_rpart_to_datatree(tree_rpart_fit)
originalTree <- Clone(causalTreeDatatree)
allFeatures <- c('H1', 'H2', 'H3', "H4", "Z")
allFeatureRanges <- list(H1 = c(lower = -Inf, upper = Inf), 
                         H2 = c(lower = -Inf, upper = Inf),
                         H3 = c(lower = -Inf, upper = Inf),
                         H4 = c(lower = -Inf, upper =Inf),
                         Z  = c(lower = 0, upper = 1))
names(allFeatureRanges) <- allFeatures
isZCategorical <- F
isHCategorical <- c(F,F,F,F)
isFeaturesCategorical <- c(isHCategorical, isZCategorical)
names(isFeaturesCategorical) <- allFeatures
  
# compute optimal Z
Z <- "Z"
varsOfInterest <- setdiff(allFeatures, Z)
generating_cohorts_by_removingZ(causalTreeDatatree, Z, allFeatures, 
                                isFeaturesCategorical, allFeatureRanges)
optimalZ_for_Xcohorts <- get_optimal_Z_for_Xcohorts(causalTreeDatatree, Z, 
                                                      varsOfInterest, originalTree,
                                                      isFeaturesCategorical,
                                                      allFeatureRanges)
# compute ATE at optimal Z
ATE_test_GCT <- ATE_at_optimalZ_continuousT_testing(testingDataFeatures, 
                                                              optimalZ_for_Xcohorts,
                                                              isFeaturesCategorical)
  
ATE_test_GCT

## Baseline1: Causal Tree-Binary (CT-B)
tree_rpart_fit_binary <- causalTree(Y1~H1+H2+H3+H4,
                                    data = trainingDataProcessed, 
                                    treatment = trainingDataProcessed$treatment,
                                    split.Rule = "CT", cv.option = "CT", 
                                    split.Honest = T, cv.Honest = T, minsize = 20)
causalTreeDatatree_binary <- convert_rpart_to_datatree(tree_rpart_fit_binary)

get_Xcohorts_from_tree_withCausalEffect <- function(treenode, varsOfInterest,
                                                      isFeaturesCategorical,
                                                      allFeatureRanges)
{
    # Obtain the Xcohorts from the tree leaves
    # Parameters:
    # -- treenode: tree node of type 'data.tree' (usually the root node)
    # -- varsOfInterest: a vector containing all features
    # -- isFeaturesCategorical: a vector indicating if the features are categorical
    # -- allFeatureRanges: a list containing the range of all features
    # Output: a list of list consisting of all Xcohorts obtained from the tree leaves
    regions <- lapply(treenode$leaves,
                      function(node) return(list(cohortRegion = get_regions_from_path(node$rpartPath,
                                                                                      varsOfInterest,
                                                                                      isFeaturesCategorical,
                                                                                      allFeatureRanges), 
                                                 causalEffect = node$causalEffect$causalEffect)))
    names(regions) <- sapply(treenode$leaves, function(node) return(node$name))
    return(regions)
}
  
Xcohorts_binary <- get_Xcohorts_from_tree_withCausalEffect(causalTreeDatatree_binary,
                                                           varsOfInterest,
                                                           isFeaturesCategorical[c(varsOfInterest)],
                                                           allFeatureRanges[c(varsOfInterest)])
  
ATE_test_CTB <- ATE_binary_CT_continuousT_testing(testingDataFeatures, 
                                                     Xcohorts_binary,
                                                     isFeaturesCategorical)
ATE_test_CTB
  
## Baseline2: Causal Tree-Multiple (CT-M)
tree_rpart_fit_binary_bins <- causalTree(Y1~H1+H2+H3+H4,
                                           data = trainingDataProcessedBins, 
                                           treatment = trainingDataProcessedBins$treatment,
                                           split.Rule = "CT", cv.option = "CT", 
                                           split.Honest = T, cv.Honest = T, minsize = 20)
  
  
causalTreeDatatree_binary_bins <- convert_rpart_to_datatree(tree_rpart_fit_binary_bins)
  
allFeatureRanges <- list(H1 = c(lower = -Inf, upper = Inf), 
                         H2 = c(lower = -Inf, upper = Inf),
                         H3 = c(lower = -Inf, upper = Inf),
                         H4 = c(lower = -Inf, upper =Inf),
                         Z = c(0.125,0.375,0.625,0.875))
names(allFeatureRanges) <- allFeatures
isZCategorical <- T
isHCategorical <- c(F,F,F,F)
isFeaturesCategorical <- c(isHCategorical, isZCategorical)

names(isFeaturesCategorical) <- allFeatures
  
Xcohorts_binary_bins <- get_Xcohorts_from_tree_withCausalEffect(causalTreeDatatree_binary_bins,
                                                                  varsOfInterest,
                                                                  isFeaturesCategorical[c(varsOfInterest)],
                                                                  allFeatureRanges[c(varsOfInterest)])
  
get_Xcohorts_from_tree_multiple_treatments_withCausalEffect <- function(treenode, varsOfInterest, Z,
                                                                          isFeaturesCategorical,
                                                                          allFeatureRanges, 
                                                                          trainingDataProcessed)
{
    
    dataList = lapply(1:nrow(trainingDataProcessed), function(i) trainingDataProcessed[i,])
    trainingDataLocation <- t(mapply(whichCohort, dataList, list(Xcohorts_binary_bins), list(isFeaturesCategorical)))
    
    get_multiple_treatment_Effects <- function(cohortId)
    {
      sample <- trainingDataProcessed[which(trainingDataLocation == cohortId),]
      Zrange <- as.numeric(get_feature_range(Z, allFeatureRanges))
      causalEffectAllZregions <- c()
      for(zlevelsId in 1:length(Zrange))
      {
        currZlevel <- Zrange[zlevelsId]
        sampleTr <- sample[which(sample$treatment==1 & sample$Z==currZlevel), "Y1"]
        sampleTe <- sample[which(sample$treatment==0 & sample$Z==currZlevel), "Y1"]
        if(length(sampleTr) > 0 && length(sampleTe) > 0)
        { causalEffect <- mean(sampleTr) - mean(sampleTe) }
        else{ causalEffect <- treenode$leaves[[cohortId]]$causalEffect$causalEffect}
        causalEffectCurrZregions <- c(Zlevel = currZlevel, causalEffect = causalEffect)
        causalEffectAllZregions <- rbind(causalEffectAllZregions, causalEffectCurrZregions)
      }
      rownames(causalEffectAllZregions) <- NULL
      return(causalEffectAllZregions)
    }
    
    numLeaves <- length(treenode$leaves)
    regions <- lapply(1:numLeaves,
                      function(nodeId) 
                      {
                        node <- treenode$leaves[[nodeId]]
                        causalEffectDF <- get_multiple_treatment_Effects(nodeId)
                        maxId <- which.max(causalEffectDF[,2])
                        return(list(cohortRegion = get_regions_from_path(node$rpartPath,
                                                                         varsOfInterest,
                                                                         isFeaturesCategorical,
                                                                         allFeatureRanges),
                                    causalEffectAllZregions = causalEffectDF,
                                    optimalCausalEffect =  causalEffectDF[maxId,2]
                        ))
                      })
    names(regions) <- sapply(treenode$leaves, function(node) return(node$name))
    
    return(regions)
}
  
Xcohorts_multiple_bins <- get_Xcohorts_from_tree_multiple_treatments_withCausalEffect(causalTreeDatatree_binary_bins,
                                                                                        varsOfInterest, Z,
                                                                                        isFeaturesCategorical,
                                                                                        allFeatureRanges, trainingDataProcessedBins)
  
  
## evaluation of binary causal tree
ATE_multiple_CT_categoricalT_testing <- function(givenH, Xcohorts_multiple, 
                                                   isFeaturesCategorical)
{
    Ntest <- dim(givenH)[1]
    numCohort <- length(Xcohorts_multiple)
    dataList = lapply(1:nrow(givenH), function(i) givenH[i,])
    findCohort <- t(mapply(whichCohort, dataList, list(Xcohorts_multiple),
                           list(isFeaturesCategorical)))
    
    get_optimal_treatment_multiple_categoricalT <- function(cohortId, Xcohorts_multiple)
    {
      causalEffectAllZLevels<- Xcohorts_multiple[[cohortId]]$causalEffectAllZregions
      optimalCausalEffect <- max(causalEffectAllZLevels[, 'causalEffect'])
      if(optimalCausalEffect <= 0) {return(0)}
      else{
        optimalZRange <- causalEffectAllZLevels[which(causalEffectAllZLevels[,'causalEffect'] == optimalCausalEffect), 'Zlevel']
        Treatment <- sample(optimalZRange, size=1)
        return(Treatment)
      }
    }
    
    treatmentForTesting <- mapply(get_optimal_treatment_multiple_categoricalT,
                                  findCohort, list(Xcohorts_multiple))
    A <- treatmentForTesting
    U.Y <- rnorm(Ntest, mean = 0, sd = uncertaintyW)
    Y <- meanFun(givenH, A) + U.Y
    return(mean(Y))
}
  
ATE_test_CTM <- ATE_multiple_CT_categoricalT_testing(testingDataFeatures, 
                                                            Xcohorts_multiple_bins,
                                                            isFeaturesCategorical)

ATE_test_CTM
################################################################################

################################################################################
## Task 2: Accuracy of estimation: MSE over training data
groundTruthEffect <- groundTruthTreatmentEffect(trainingDataProcessed)
id <- which(trainingDataProcessed$A != 0)
groundTruthEffectTreatments <- groundTruthEffect[id]

## GCT
trainingDataProcessedBack <- data.frame(trainingDataProcessed[id,c('H1', 'H2', 'H3', 'H4', 'treatment', 'A')])
colnames(trainingDataProcessedBack) <- c('H1', 'H2', 'H3', 'H4', 'treatment', 'Z')
estimatedTreatmentEffect <- predict(tree_rpart_fit, trainingDataProcessedBack)
MSE_GCT <- mean((estimatedTreatmentEffect-groundTruthEffectTreatments)**2)
MSE_GCT
  
## Baseline1: CTB
estimatedTreatmentEffectBinary <- predict(tree_rpart_fit_binary, trainingDataProcessed[id,])
MSE_CTB <- mean((estimatedTreatmentEffectBinary-groundTruthEffectTreatments)**2)
MSE_CTB
  
## Baseline2: CTM 
numCohort <- length(Xcohorts_multiple_bins)
dataList = lapply(1:nrow(trainingDataProcessedBins[id,]),
                    function(i) trainingDataProcessedBins[i,])
findCohort <- t(mapply(whichCohort, dataList, list(Xcohorts_multiple_bins),
                         list(isFeaturesCategorical)))
estimatedTreatmentEffectMultiple <- 
    sapply(id,
           function(sampleId) {
             currCohort <- Xcohorts_multiple_bins[[findCohort[sampleId]]]$causalEffectAllZregions;
             treatmentReceived <- trainingDataProcessedBins[sampleId,'A']
             est <- currCohort[which(currCohort[,'Zlevel'] == treatmentReceived), 'causalEffect']
             return(est)
           }
    )
MSE_CTM <- mean((estimatedTreatmentEffectMultiple-groundTruthEffectTreatments)**2)
MSE_CTM 
################################################################################

################################################################################ 
## Output (continuous treatments)
output_ATE_test <- c(ATE_test_GCT, ATE_test_CTB, ATE_test_CTM)
output_ATE_test

output_MSE <- c(MSE_GCT, MSE_CTB, MSE_CTM)
output_MSE
################################################################################ 


