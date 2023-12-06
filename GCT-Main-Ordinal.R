################################################################################
##               Generalized Causal Tree (GCT) for Uplift Modeling            ## 
##                           (Ordinal treatments)                             ##
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
source('GCT-Utils-Simulation-Ordinal.R')

################################################################################
## Generate simulated data 
seed = 2
repId = 1
w1 = rnorm(1)
w2 = rnorm(1)
uncertaintyW <<- 1
suppressMessages(suppressWarnings(
  {
    trainingData <- generate_ordinalT_training(seed = seed+repId, size = 500, 
                                uncertaintyW=uncertaintyW, w1, w2)
    trainingDataProcessed <- process_nonbinary_treatment(trainingData) 
    testingData <- generate_ordinalT_training(seed = seed+repId+1, size = 500,
                               uncertaintyW=uncertaintyW, w1, w2)
    testingDataFeatures <- testingData[, c('ID', 'H1', 'H2', 'H3', 'H4')]
  }))
  
trainingDataExport <- trainingData[, c('H1', 'H2', 'A', 'Y1')]
testingDataExport <- testingData[, c('H1', 'H2')]
################################################################################
  
################################################################################  
## Task 1: Maximizing the average outcome with optimal treatment allocations
## over the testing data
  
## Method: GCT
tree_rpart_fit <- causalTree(Y1~H1+H2+H3+H4+factor(Z),
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
                         Z = c(1,2,3,4,5,6))
names(allFeatureRanges) <- allFeatures
isZCategorical <- T
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
ATE_test_GCT <- ATE_at_optimalZ_categoricalT_testing(testingDataFeatures, 
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
  
ATE_test_CTB <- ATE_binary_CT_categoricalT_testing(testingDataFeatures, 
                                                       Xcohorts_binary,
                                                       isFeaturesCategorical)
ATE_test_CTB

## Baseline2: Causal Tree-Multiple (CT-M)
get_Xcohorts_from_tree_multiple_treatments_withCausalEffect <- function(treenode, varsOfInterest, Z,
                                                                          isFeaturesCategorical,
                                                                          allFeatureRanges, 
                                                                          trainingDataProcessed)
  {
    
    dataList = lapply(1:nrow(trainingDataProcessed), function(i) trainingDataProcessed[i,])
    trainingDataLocation <- t(mapply(whichCohort, dataList, list( Xcohorts_binary), list(isFeaturesCategorical)))
    
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
  
Xcohorts_multiple <- get_Xcohorts_from_tree_multiple_treatments_withCausalEffect(causalTreeDatatree_binary,
                                                                                    varsOfInterest, Z,
                                                                                    isFeaturesCategorical,
                                                                                    allFeatureRanges, trainingDataProcessed)
  
ATE_test_CTM <- ATE_multiple_CT_categoricalT_testing(testingDataFeatures, 
                                                        Xcohorts_multiple,
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
estimatedTreatmentEffectBinary <- predict(tree_rpart_fit_binary, 
                                            trainingDataProcessed[id,],
                                            type = "vector")
MSE_CTB <- mean((estimatedTreatmentEffectBinary-groundTruthEffectTreatments)**2)
MSE_CTB
  
## Baseline2: CTM
numCohort <- length(Xcohorts_multiple)
dataList = lapply(1:nrow(trainingDataProcessed[id,]),
                    function(i) trainingDataProcessed[i,])
findCohort <- t(mapply(whichCohort, dataList, list(Xcohorts_multiple),
                         list(isFeaturesCategorical)))
estimatedTreatmentEffectMultiple <- 
    sapply(id,
           function(sampleId) {
             currCohort <- Xcohorts_multiple[[findCohort[sampleId]]]$causalEffectAllZregions;
             treatmentReceived <- trainingDataProcessed[sampleId,'A']
             est <- currCohort[which(currCohort[,'Zlevel'] == treatmentReceived), 'causalEffect']
             return(est)
           }
    )
  MSE_CTM <- mean((estimatedTreatmentEffectMultiple-groundTruthEffectTreatments)**2)
  MSE_CTM
  
}
################################################################################

################################################################################ 
## Output (categorical treatments)
output_ATE_test <- c(ATE_test_GCT, ATE_test_CTB, ATE_test_CTM)
output_ATE_test

output_MSE <- c(MSE_GCT, MSE_CTB, MSE_CTM)
output_MSE
################################################################################ 
