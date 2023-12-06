# GCT-Utils-OptimalCohorts.R
traverse_levelorder <- function(treenode, varsOfInterest,
                                isFeaturesCategorical, allFeatureRanges)
{
  # Level-order traversal on a binary tree of type 'data.tree'
  # to remove empty branches
  # Parameters:
  # -- treenode: tree node of type 'data.tree' (usually the root node)
  # -- varsOfInterest: a vector containing all features
  # -- isFeaturesCategorical: a vector indicating if the features are categorical
  # -- allFeatureRanges: a list containing the range of all features
  # Output: NULL
  for(i in 1:treenode$height)
  {
    traverse_levelorder_given_level(treenode, i, varsOfInterest,
                                    isFeaturesCategorical, allFeatureRanges)
  }
}

traverse_levelorder_given_level <- function(treenode, level, varsOfInterest,
                                            isFeaturesCategorical, allFeatureRanges)
{
  # Level-order traversal on a binary tree of type 'data.tree'
  # for given level to remove empty branches
  # Parameters:
  # -- treenode: tree node of type 'data.tree'
  # -- level: level of interest
  # -- varsOfInterest: a vector containing all features
  # -- isFeaturesCategorical: a vector indicating if the features are categorical
  # -- allFeatureRanges: a list containing the range of all features
  # Output: NULL
  if(!is.null(treenode))
  {
    if(level == 1){
      #print(paste("Current node:", treenode$name))
      currentRegion <- get_regions_from_path(treenode$rpartPath,
                                             varsOfInterest,
                                             isFeaturesCategorical,
                                             allFeatureRanges)
      if(!isRegionValid(currentRegion))
      {
        # if the current region is not valid
        if(treenode$position==1)
        {
          parName <- treenode$parent$name
          if(any(grep(">=", parName))) {newparName <- sub(">=", "<", parName)}
          else if(any(grep(">", parName)))
          {newparName <- sub(">", "<=", parName)}
          else if(any(grep("<=", parName)))
          {newparName <- sub("<=", ">", parName)}
          else if(any(grep("<", parName)))
          {newparName <- sub("<", ">=", parName)}
          treenode$parent$name <- newparName
        }
        treenode$parent$RemoveChild(treenode$name)
      }
    }
    else{
      if(length(treenode$children) == 1)
      { traverse_levelorder_given_level(treenode$children[[1]],
                                        level-1, varsOfInterest,
                                        isFeaturesCategorical,
                                        allFeatureRanges)}
      else if(length(treenode$children) == 2)
      {
        lchild <- treenode$children[[1]]
        rchild <- treenode$children[[2]]
        traverse_levelorder_given_level(lchild, level-1, varsOfInterest,
                                        isFeaturesCategorical, allFeatureRanges)
        traverse_levelorder_given_level(rchild, level-1, varsOfInterest,
                                        isFeaturesCategorical, allFeatureRanges)
      }
      else if(length(treenode$children) > 2)
      {
        print("This node has more than 2 children.")
        childrens <- treenode$children
        for(i in 1:length(childrens))
        {traverse_levelorder_given_level(childrens[[i]],level-1,varsOfInterest,
                                         isFeaturesCategorical,
                                         allFeatureRanges)}
      }
    }
  }
}

isIntervalValid <- function(interval)
{
  # Check whether an interval is valid or not
  # Parameters:
  # -- interval: c(lower,upper)
  # Output: TRUE/FALSE
  ## Old version
  #if(length(interval) == 0) return(FALSE)
  #else if (interval[1] > interval[2]) return(FALSE)
  #else return(TRUE)
  
  if(length(interval) == 0) return(FALSE)
  else if (length(interval) == 2 && (interval[1] > interval[2])) return(FALSE)
  else return(TRUE)
}


isRegionValid <- function(region)
{
  # Check whether a region is valid or not
  # Parameters:
  # -- region: a list, each element is a interval or a set.
  # Output: TRUE/FALSE
  numFeature <- length(region)
  isRegionValidVec <- sapply(1:numFeature,
                             function(i) isIntervalValid(region[[i]]))
  if(any(!isRegionValidVec)) return(FALSE)
  else return(TRUE)
}

removeZ <- function(treenode, Z)
{
  # Remove nodes associated with the variable name Z
  # Parameters:
  # -- treenode: tree node of type 'data.tree' (usually the root node)
  # -- Z: the variable to be removed
  # Output: NULL
  swapChildren <- function(treenode)
    ## swap the position of two children
  {
    if (length(treenode$children) == 2)
    {
      childl <- treenode$children[[1]]
      childr <- treenode$children[[2]]
      childClone <- Clone(childl)
      treenode$RemoveChild(childl$name)
      treenode$AddChildNode(childClone)
    }
  }
  ZinPath <- grep(Z, treenode$rpartPath)
  if(any(ZinPath)) { treenode$rpartPath <- treenode$rpartPath[-tail(ZinPath,1)]}
  ## remove Z from the recorded path information
  if(!treenode$isLeaf)
  {
    for(j in 1:length(treenode$children))
    {
      ## post-order traversal using recursion
      removeZ(treenode$children[[j]],Z)
    }
    # print(paste("Current Node:", treenode$name))
    if(treenode$nodeName == Z || treenode$nodeName == paste0('factor(', Z, ')'))
    {
      print(paste("Remove node: ", treenode$name))
      children <- treenode$children
      numChildren <- length(children)
      if(numChildren == 1)
      {
        ## the node has only one child
        treenode$parent$AddChildNode(children[[1]])
        treenode$parent$RemoveChild(treenode$name)
      }
      else{
        ## the node has two children
        childIdWrtItsPar <- treenode$position
        ## record the information on whether the current Z node is the first
        ## or second child w.r.t. its parent
        isLeafVec <- sapply(1:2, function(i) return(children[[i]]$isLeaf))
        if(all(isLeafVec))
          ## if both children are leaf node
        {
          ## union
          mergeLeafNodeName <- paste0("( ", children[[1]]$name, "\n U ",
                                      children[[2]]$name," )")
          treenode$parent$AddChild(mergeLeafNodeName,
                                   rpartPath = treenode$rpartPath,
                                   nodeName = 'leaf',
                                   causalEffect = rbind(
                                     children[[1]]$causalEffect,
                                     children[[2]]$causalEffect))
          treenodePar <- treenode$parent
          treenode$parent$RemoveChild(treenode$name)
          if(childIdWrtItsPar == 1)
          { swapChildren(treenodePar) }
          print("successful removal")
        }
        else if(any(isLeafVec))
          ## one child is leaf node but the other is not
        {
          idNotLeaf <- which(!isLeafVec)
          idLeaf <- which(isLeafVec)
          removeNodeName <- children[[idLeaf]]$name
          removeNodeCausalEffect <-  children[[idLeaf]]$causalEffect
          treenodePar <- treenode$parent
          mergeLeafNodeName <-
            lapply(children[[idNotLeaf]]$leaves,
                   function(x) {
                     x$name <-  paste(x$name," \n + ",removeNodeName);
                     x$causalEffect <- rbind(x$causalEffect,
                                             removeNodeCausalEffect)})
          treenode$parent$AddChildNode(children[[idNotLeaf]])
          treenode$parent$RemoveChild(treenode$name)
          if(childIdWrtItsPar == 1)
          { swapChildren(treenodePar) }
          print("successful removal")
        }
        else
          ## both children are not leaf nodes
        {
          ## merge two subtrees
          leafCountSubtree1 <- children[[1]]$leafCount
          leafCountSubtree2 <- children[[2]]$leafCount
          ## add tree1 as a child of treenode$parent
          ## add tree2 to each leaf node of tree1
          treenode$parent$AddChildNode(children[[1]])
          treenodePar <- treenode$parent
          treenode$parent$RemoveChild(treenode$name)
          if(childIdWrtItsPar == 1)
          { swapChildren(treenodePar) }
          leavesSubtree1 <- children[[1]]$leaves
          leavesNames <- sapply(1:leafCountSubtree1, function(i)
            leavesSubtree1[[i]]$name)
          for(i in 1:leafCountSubtree1)
          {
            leaf <- FindNode(children[[1]]$parent, leavesNames[i])
            leafPar <- leaf$parent
            leafPos <- leaf$position
            cloneTree2 <- Clone(children[[2]])
            mergeLeafNodeName <- lapply(cloneTree2$leaves,
                                        function(x)
                                        { x$name <- paste(leaf$name,
                                                          "\n + " ,x$name);
                                        x$causalEffect<-rbind(leaf$causalEffect,
                                                              x$causalEffect)})
            ## change the $rpartPath information of cloneTree2
            mergerpartPath <- Traverse(cloneTree2, traversal = "pre-order",
                                       filterFun = function(x)
                                       {x$rpartPath <- c(leaf$rpartPath,
                                                         setdiff(x$rpartPath,
                                                                 cloneTree2$rpartPath));
                                       return(TRUE)})
            cloneTree2$name <- paste0(cloneTree2$name, ":clone", i)
            leafPar$AddChildNode(cloneTree2)
            ## We cannot rename the node here, because in data.tree object,
            ## the nodes are uniquely identified by their names.
            ## It is not allowed to have two nodes that shares exactly the
            ## same node in a data.tree.
            #renameNode <- FindNode(leafPar, cloneTree2$name)
            #renameNode$name <- children[[2]]$name
            leafPar$RemoveChild(leavesNames[i])
            if(leafPos == 1) { swapChildren(leafPar) }
          }
          print("successful removal")
        }
      }
    }
  }
}

generating_cohorts_by_removingZ <- function(treenode, Z, allFeatures,
                                            isFeaturesCategorical,
                                            allFeatureRanges)
{
  # Generate cohorts after removing Z and removing empty branches
  # Parameters:
  # -- treenode: tree node of type 'data.tree' (usually the root node)
  # -- Z: the variable to be removed
  # -- varsOfInterest: a vector containing all features
  # -- isFeaturesCategorical: a vector indicating if the features are categorical
  # -- allFeatureRanges: a list containing the range of all features
  # Output: NULL
  removeZ(treenode,Z)
  varsOfInterest <- setdiff(allFeatures, Z)
  traverse_levelorder(treenode, varsOfInterest,
                      isFeaturesCategorical, allFeatureRanges)
}

get_feature_range <- function(feature, allFeatureRanges)
{
  # Obtain the range for a given feature
  # Parameters:
  # -- feature: feature name;
  # -- allFeatureRanges: a list containing the range of all features;
  # Output: the range for a given feature;
  #         for a continuous variable, the output is an interval
  #         for a categorical variable, the output is a vector
  featureId <- which(names(allFeatureRanges) == feature)
  if(length(featureId) > 0) {return(allFeatureRanges[[featureId]])}
  else {message("The range of this feature was not specified.")}
}

get_feature_type <- function(feature, isFeaturesCategorical)
{
  # Obtain the type (categorical or not) for a given feature
  # Parameters:
  # -- feature: feature name;
  # -- isFeaturesCategorical: a vector indicating if the features are categorical;
  # Output: TRUE/FALSE
  featureId <- which(names(isFeaturesCategorical) == feature)
  if(length(featureId) > 0) {return(isFeaturesCategorical[[featureId]])}
  else {message("The type of this feature was not specified.")}
}

convert_text_to_interval <- function(textList, feature,
                                     isFeatureCategorical,
                                     featureRange)
{
  # Convert mathematical expressions into an interval
  # Parameters:
  # -- text: mathematical expression of inequality (type: a vector of string)
  # -- feature: name of a give feature  (type: string)
  # -- isFeatureCategorical: a logic value indicating if
  #                          the feature of interest is categorical
  # -- featureRange: specify the range of the feature of interest;
  #                  for a continuous variable, the input is an interval
  #                  for a categorical variable, the input is a vector
  #                  consisting of all possible categories.
  # Output: for a continuous variable, the output is an interval;
  #         for a categorical variable, the output is a vector (representing a set)
  # eg1: convert_text_to_interval('x1>= -0.9932', 'x1')
  # eg2: convert_text_to_interval('x1< 0.9932', 'x1')
  # eg3: convert_text_to_interval('x1 = 0.9932', 'x1')
  # eg4: convert_text_to_interval(c("x1< -0.06912", "x1>=-0.9932"), 'x1')
  # eg5: convert_text_to_interval(c("x1<2","x1>=-9","x1>=-8","x1<=-1"),'x1')
  # eg6: convert_text_to_interval(c("factor(Z)=1,2,3", "factor(Z)=1,2"), 'Z',
  #                              isFeatureCategorical = TRUE,
  #                              featureRange = c(1,2,3,4))
  if(!isFeatureCategorical){
    interval <- featureRange
    ## Initialize the interval to be (-Inf, Inf)
    names(interval) <- c('lower', 'upper')
    textListLen <- length(textList)
    if(textListLen == 0){
      return(inverval)
    }
    else{
      for(textId in 1:textListLen)
      {
        text <- textList[textId]
        textRmName <- str_remove(text, feature)
        ## remove characters corresponding to the 'feature' from the 'text' string
        num <- parse_number(textRmName)
        ## extract the number from the 'text' string
        if(any(grep(">", text)))
        {
          ## if the expression is of the form of 'X > num',
          ## then the interval is (a, Inf)
          interval[1] <- max(num, interval[1])
        }
        else if (any(grep("<", text)))
        {
          ## if the expression is of the form of 'X < num',
          ## then the interval is (-Inf, num)
          interval[2] <- min(num, interval[2])
        }
        else if (any(grep("=", text)))
        {
          ## if the expression is of the form of 'X = num',
          ## then the interval is (num, num)
          interval[1] <- num
          interval[2] <- num
        }
        else
        {
          return("Error in converting text to interval!")
        }
      }
      return(interval)
    }
  }
  else{
    set <- featureRange
    textListLen <- length(textList)
    if(textListLen == 0){
      return(set)
    }
    else{
      for(textId in 1:textListLen)
      {
        text <- textList[textId]
        textRmName <- str_remove(text, fixed(paste0('factor(', feature, ')=')))
        setElements <- strsplit(textRmName, ',')[[1]]
        set <- intersect(set, setElements)
      }
      return(set)
    }
  }
}

get_regions_from_path_given_feature <- function(splitPath, feature,
                                                isFeatureCategorical,
                                                featureRange)
{
  # obtain regions of variables from one split paths for given feature
  # Parameters:
  # -- splitPath: list of splits
  # -- feature: variable name
  # -- isFeatureCategorical: a logic value indicating if
  #                          the feature of interest is categorical
  # -- featureRange: specify the range of the feature of interest;
  #                  for a continuous variable, the input is an inverval
  #                  for a categorical variable, the input is a vector
  #                  consisting of all possible categories.
  # Output: for a continuous variable, the output is an interval;
  #         for a categorical variable, the output is a vector (representing a set)
  pathLen <- length(splitPath)
  if(length(grep(feature, splitPath)) != pathLen)
    return("The input contains unmatched features!")
  region_given_feature <- convert_text_to_interval(splitPath, feature,
                                                   isFeatureCategorical,
                                                   featureRange)
  return(region_given_feature)
}

get_regions_from_path <- function(splitPath, varsOfInterest,
                                  isFeaturesCategorical, allFeatureRanges)
{
  # obtain regions of variables from split paths
  # Parameters:
  # -- splitPath: vector of splits
  # -- varsOfInterest: all features
  # -- isFeaturesCategorical: a vector indicating if the features are categorical
  # -- allFeatureRanges: a list containing the range of all features
  # Output: a list consisting of regions for all variables of interest for a given path
  regionsSingleNode = c()
  for(varID in 1:length(varsOfInterest))
  {
    feature <- varsOfInterest[varID]
    varInPathId <- grep(feature, splitPath)
    if(any(varInPathId))
    {
      region_given_feature <-
        get_regions_from_path_given_feature(
          splitPath[varInPathId], feature,
          get_feature_type(feature, isFeaturesCategorical),
          get_feature_range(feature, allFeatureRanges))
      regionsSingleNode = c(regionsSingleNode, list(region_given_feature))
    }
    else
    {
      region_given_feature <- get_feature_range(feature, allFeatureRanges)
      regionsSingleNode <- c(regionsSingleNode, list(region_given_feature))
    }
  }
  names(regionsSingleNode) <- varsOfInterest
  return(regionsSingleNode)
}

get_Xcohorts_from_tree <- function(treenode, varsOfInterest,
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
                    function(node) return(get_regions_from_path(node$rpartPath,
                                                                varsOfInterest,
                                                                isFeaturesCategorical,
                                                                allFeatureRanges)))
  names(regions) <- sapply(treenode$leaves, function(node) return(node$name))
  return(regions)
}

get_optimal_Z_for_single_Xcohort <-
  function(treenode, Z, varsOfInterest, originalTree,
           isFeaturesCategorical, allFeatureRanges)
  {
    # Obtain the optimal treatment level Z based on the estimated causal effects
    # for each X cohort after removing Z nodes from the original causal tree
    # Parameters:
    # -- treenode: tree node of type 'data.tree' (usually the root node)
    # -- Z: the variable to be removed
    # -- varsOfInterest: variable list
    # -- originalTree: the original tree before removing Z nodes
    # -- isFeaturesCategorical: a vector indicating if the features are categorical
    # -- allFeatureRanges: a list containing the range of all features
    # Output: a list consisting of
    # -- $cohortName: the name of the current region
    #                (to keep track of the source leaf nodes in the original tree)
    # -- $cohortRegion: regions w.r.t. X variables
    # -- $causalEffectAllZregions: the causal effects extracted from the source leaf nodes
    # -- $optimalCausalEffect: the optimal causal effects for this X cohort 
    # -- $optimalZregion: the Z levels to which the optimal causal effects corresponds
    causalEffectDF <- treenode$causalEffect
    maxId <- which.max(causalEffectDF[,2])
    maxLeaf <- causalEffectDF[maxId, 1]
    matchLeafId <- which(sapply(originalTree$leaves, function(x) x$name)==maxLeaf)
    return(list(
      cohortName = treenode$name,
      cohortRegion = get_regions_from_path(
        treenode$rpartPath, varsOfInterest,
        isFeaturesCategorical, allFeatureRanges),
      causalEffectAllZregions =  causalEffectDF,
      optimalCausalEffect = causalEffectDF[maxId,2],
      optimalZregion = get_regions_from_path(
        originalTree$leaves[[matchLeafId]]$rpartPath, Z,
        isFeaturesCategorical, allFeatureRanges)))
  }

get_optimal_Z_for_Xcohorts <- function(treenode, Z, varsOfInterest,
                                       originalTree, isFeaturesCategorical,
                                       allFeatureRanges)
{
  # Obtain the optimal treatment level Z based on the estimated causal effects
  # for all X cohorts after removing Z nodes from the original causal tree
  # Parameters:
  # -- treenode: tree node of type 'data.tree' (usually the root node)
  # -- Z: the variable to be removed
  # -- varsOfInterest: variable list
  # -- originalTree: the original tree before removing Z nodes
  # -- isFeaturesCategorical: a vector indicating if the features are categorical
  # -- allFeatureRanges: a list containing the range of all features
  # Output: a list of list consisting of outputs from get_optimal_Z_for_single_Xcohort()
  #         for all Xcohorts in the tree
  Xcohorts <- treenode$leaves
  optimalZlist <- c()
  for(i in 1:length(Xcohorts))
  {
    optimalZlist = c(optimalZlist,
                     list(get_optimal_Z_for_single_Xcohort(
                       Xcohorts[[i]], Z, varsOfInterest, originalTree,
                       isFeaturesCategorical, allFeatureRanges)))
  }
  return(optimalZlist)
}

get_minimal_Z_from_tree <- function(Z, OriginalTree, isFeaturesCategorical,
                                    allFeatureRanges)
{
  # Obtain the minimal partitions of treatment levels Z from the original tree
  # Parameters:
  # -- originalTree: the original tree before removing Z nodes
  # -- Z: the variable to be removed
  # -- isFeaturesCategorical: a vector indicating if the features are categorical
  # -- allFeatureRanges: a list containing the range of all features
  # Output: for a continuous variable, the output is an matrix
  #         (each row is an interval indicating range for one treatment partition);
  #         for a categorical variable, the output is a vector;
  Zrange <- get_feature_range(Z, allFeatureRanges)
  isZcategorical <- get_feature_type(Z, isFeaturesCategorical)
  if(isZcategorical) # Z is a categorical variable
  {
    return(Zrange)
  }
  else # Z is a continuous variable
  {
    regions <- sapply(originalTree$leaves,
                      function(node) return(get_regions_from_path(node$rpartPath, Z,
                                                                  isFeaturesCategorical,
                                                                  allFeatureRanges)))
    regions <- unique(regions)
    if(length(regions) == 1) return(regions)
    cutPoints <- sort(unique(unlist(lapply(regions, function(x) as.vector(x)))))
    minimalZregions <- matrix(c(head(cutPoints,-1),
                                head(shift(cutPoints, type = "lead"),-1)), ncol=2)
    colnames(minimalZregions) <- c('lower', 'upper')
    return(minimalZregions)
  }
}

is_Zpartition_in_Zregion <- function(Zpartition, Zregion, isZcategorical)
{
  # Determine whether a given Zpartition is in the given Zregion
  # Parameters:
  # -- isZcategorical: TRUE/FALSE indicating whether Z is categorical
  # -- If Z is a categorical variable:
  # ---- Zpartition: a charactor/number specifying the Z level;
  # ---- Zregion: a vector (respresenting a set)
  # -- If Z is a continuous variable:
  # ---- Zpartition: an interval
  # ---- Zregion: an interval
  # Output: TRUE/FALSE
  if(isZcategorical)
  {
    return(Zpartition %in% Zregion)
  }
  else{
    return( Zpartition[1] >= Zregion[1] && Zpartition[2] <= Zregion[2])
  }
}

get_effects_for_Xcohorts_allZ <- function(treenode, Z, varsOfInterest,
                                          originalTree, isFeaturesCategorical,
                                          allFeatureRanges)
{
  # Obtain the estimated causal effects treatment of each level Z
  # for all X cohorts after removing Z nodes
  # Parameters:
  # -- treenode: tree node of type 'data.tree' (usually the root node)
  # -- Z: the variable to be removed
  # -- varsOfInterest: variable list
  # -- originalTree: the original tree before removing Z nodes
  # -- isFeaturesCategorical: a vector indicating if the features are categorical
  # -- allFeatureRanges: a list containing the range of all features
  # Output: a list of list consisting of
  # -- $cohortName: the name of the current region
  # -- $cohortRegion: regions w.r.t. X variables
  # -- $causalEffectForEachZregion: the estimated causal effects treatment of each
  #                                 treatment level Z within this X cohort
  isZcategorical <- get_feature_type(Z, isFeaturesCategorical)
  minimalZregions <- get_minimal_Z_from_tree(Z, originalTree,
                                             isFeaturesCategorical, allFeatureRanges)
  numMinimalZregions <- ifelse(isZcategorical, length(minimalZregions),
                               dim(minimalZregions)[1])
  Xcohorts <- treenode$leaves
  effects_for_Xcohorts_allZlist <- c()
  for(i in 1:length(Xcohorts))
  {
    currCohort <- Xcohorts[[i]]
    causalEffectDF <- currCohort$causalEffect
    causalEffectGivenZ <- rep(0, numMinimalZregions)
    effectInEachRegion <- data.frame(minimalZregions, causalEffectGivenZ)
    for(leafId in 1:dim(causalEffectDF)[1])
    {
      matchLeafId <- which(sapply(originalTree$leaves,
                                  function(x) x$name)==causalEffectDF[leafId,1])
      matchLeafZ <- get_regions_from_path(originalTree$leaves[[matchLeafId]]$rpartPath, Z,
                                          isFeaturesCategorical, allFeatureRanges)$Z
      if(isZcategorical)
      {
        minimalMatchIndicator <- sapply(1:numMinimalZregions, function(j)
          is_Zpartition_in_Zregion(minimalZregions[j], matchLeafZ, isZcategorical))
      }
      else{
        minimalMatchIndicator <- sapply(1:numMinimalZregions, function(j)
          is_Zpartition_in_Zregion(minimalZregions[j, ], matchLeafZ, isZcategorical))
      }
      effectInEachRegion$causalEffectGivenZ[minimalMatchIndicator] <-
        causalEffectDF$causalEffect[leafId]
      effects_for_Xcohorts_allZ <- list(cohortName = currCohort$name,
                                        cohortRegion = get_regions_from_path(
                                          currCohort$rpartPath, varsOfInterest,
                                          isFeaturesCategorical, allFeatureRanges),
                                        causalEffectForEachZregion = effectInEachRegion)
    }
    effects_for_Xcohorts_allZlist <- c(effects_for_Xcohorts_allZlist,
                                       list(effects_for_Xcohorts_allZ))
  }
  return(effects_for_Xcohorts_allZlist)
}
