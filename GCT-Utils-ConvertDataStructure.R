# GCT-Utils-ConvertDataStructure.R
rpart_id_to_var_name <- function(rpart_id, rpart_tree){
  # Obtain the variable name associated with each node
  # Parameters:
  # -- rpart_id: number (int)
  # -- rpart_tree: a tree of type 'rpart'
  # Output: variable name
  row_id = which(rownames(rpart_tree$frame) == rpart_id)
  if(length(row_id) == 0)
  {
    print(paste("The Node", rpart_id, "is not in this tree!"))
    return(character(0))
  }
  else return(rpart_tree$frame$var[row_id])
}

traverse_preorder <- function(treenode, rpart_tree){
  # Pre-order traversal on a binary tree of type 'data.tree' to
  # recover the node path ($nodeName) and node variable ($rpartPath)
  # information from the 'rpart' tree.
  # Parameters:
  # -- treenode: tree node of type 'data.tree' (usually the root node)
  # -- rpart_tree: a tree of type 'rpart'
  # Output: NULL
  if(treenode$isLeaf){
    treenode$rpartPath <- path.rpart(rpart_tree, treenode$rpart.id,
                                     pretty = 0, print.it = F)[[1]]
    treenode$nodeName <- rpart_id_to_var_name(treenode$rpart.id, rpart_tree)
    return('End')
  } else {
    treenode$rpartPath <- path.rpart(rpart_tree, treenode$rpart.id,
                                     pretty = 0, print.it = F)[[1]]
    treenode$nodeName <- rpart_id_to_var_name(treenode$rpart.id, rpart_tree)
  }
  children = treenode$children
  traverse_preorder(children[[1]], rpart_tree)
  traverse_preorder(children[[2]], rpart_tree)
}

causal_effects_from_causalTree <- function(treenode)
{
  # Recover the information of the estimated causal effects ($causalEffect)
  # for leaf nodes, and rename the leaf nodes using letters.
  # Parameters:
  # -- treenode: tree node of type 'data.tree' (usually the root node)
  # Output: NULL
  leafNodes <- treenode$leaves
  for(i in 1:length(leafNodes))
  {
    leaf <- leafNodes[[i]]
    causalEffectVal <- as.numeric(leaf$name)
    if(is.na(causalEffectVal))
    {
      print("The causal effects are missing.")
      leaf$causalEffect <- 0
    }
    else
    { leaf$causalEffect <- data.frame(leafname = paste("Leaf:", i),
                                      causalEffect = causalEffectVal) }
    leaf$name <- paste("Leaf:",i)
  }
}

convert_rpart_to_datatree <- function(rpart_tree)
{
  # Convert a causal tree of type 'rpart' to type 'data.tree'
  # Parameters:
  # -- rpart_tree: a tree of type 'rpart'
  # Output: the root node of a causal tree of type 'data.tree'
  digits <- getOption("digits") - 1
  use.n  <- FALSE
  frame       <- rpart_tree$frame
  ylevels     <- attr(rpart_tree, "ylevels")
  nodes       <- as.numeric(rownames(frame))
  leaves      <- frame$var == "<leaf>"
  leaf_labels <- rpart_tree$functions$text(
    yval = if(is.null(frame$yval2)) frame$yval[leaves] else frame$yval2[leaves,],
    dev    = frame$dev[leaves],
    wt     = frame$wt[leaves],
    ylevel = ylevels,
    digits = digits,
    n      = frame$n[leaves],
    use.n  = use.n)
  nonLeaf_labels <- labels(rpart_tree)[which(!leaves) + 1L]
  while(any(duplicated(nonLeaf_labels)))
  {
    dupId <- which(duplicated(nonLeaf_labels))[1]
    nonLeaf_labels[dupId] <- paste0(nonLeaf_labels[dupId], " clone")
  }
  node_labels <- setNames(c(nonLeaf_labels, leaf_labels),
                          c(nodes[!leaves], nodes[leaves]))
  network_df  <- data.frame(
    from = node_labels[as.character(floor(nodes[-1L] / 2L))],
    to = node_labels[as.character(nodes[-1L])],
    rpart.id = nodes[-1L])
  datatree <- FromDataFrameNetwork(network_df)
  datatree$rpart.id <- nodes[1L]
  if(sum(rpart_tree$frame$var == '<leaf>') != datatree$leafCount)
  { message("error in converting the rpart tree to a data.tree tree")}
  rootNode = Node$new('RootHead', nodeName = 'rootNode')
  rootNode$AddChildNode(datatree)
  traverse_preorder(datatree, rpart_tree)
  causal_effects_from_causalTree(datatree)
  return(rootNode)
}
