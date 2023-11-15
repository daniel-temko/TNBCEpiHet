library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(cowplot)

#' Returns ids indexing 'size' random positions with 
#' each label, or all positions indexing the label, if fewer than
#' size
BalancedSample <- function(labels, size, seed = 1, replace=FALSE){
  types <- unique(labels)
  ids <- c()
  set.seed(seed)
  for(i in 1:length(types)){
    typeIds <- which(labels == types[i])
    if(length(typeIds) < size) { 
      ids <- c(ids, typeIds)
    } else {
      ids <- c(ids, sample(typeIds, size, replace = replace))
    }
  }
  return(ids)
}

#' Returns ids indexing n random positions for each label, where
#' n is given by scheme. Throws an error if the number of requested
#' items is greater than the number of possible positions when replace
#' is false
StratifiedSample <- function(labels, scheme, seed, replace=FALSE){
  types <- names(scheme)
  ids <- c()
  set.seed(seed)
  for(i in 1:length(scheme)){
    type_ids <- which(labels == types[i])
    ids <- c(ids, sample(type_ids, scheme[i], replace = replace))
  }
  return(ids)
}
