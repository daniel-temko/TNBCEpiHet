library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

#' Load RG-set information to minfi from file
read_minfi <- function(path, strip_from_names = "Bojana_"){
  tgts <- read.metharray.sheet(path)
  message("Loading data...")
  rg_set <- read.metharray.exp(targets = tgts)
  
  if(!is.na(strip_from_names)){
    message("Correcting sample names...")
    samp_names <- sapply(tgts$Sample_Name, function(x) gsub(strip_from_names, "", x))
  } else {
    samp_names <- tgts$Sample_Name
  }
  sampleNames(rg_set) <- samp_names
  
  message("Retrieving detection P-values...")
  det_p <- detectionP(rg_set)
  return(list(rg_set = rg_set, det_p = det_p))
}

#' Aggregate M-values or Beta-values by averaging to probes
#' annotated to the same region
aggregate_probes <- function(map, probe_vals){
  map <- map[which(map$locus %in% rownames(probe_vals)),]
  
  locus_vals <- probe_vals[as.character(map$locus),]
  agg_df <- aggregate(locus_vals, by = list(map$region.id), mean)
  rownames(agg_df) <- agg_df$Group.1 
  agg_df$Group.1 <- NULL
  return(agg_df)
}

#' Subset M-values or Beta-values to probes annotated to a
#' region type
cut_probes <- function(ann, probe_vals){
  ann <- ann[which(rownames(ann) %in% rownames(probe_vals)),]
  cut_df <- probe_vals[rownames(ann), ]
}

#' Aggregate M-value and Beta values using mean within regions of interest. And take cuts for
#' probe-sets of interest
aggregate_values <- function(m_vals, beta_vals, region_map_list, region_ann_list, locus_ann_list, prefix){
  return_list <- list()
  
  # Aggregate methylation values to features
  for(i in 1:length(region_map_list)){
    type <- names(region_map_list)[i]
    map <- region_map_list[[i]]
    ann <- region_ann_list[[i]]
    message("Aggregating ", type)
    
    # Aggregate
    m_agg_df <- aggregate_probes(map, m_vals)
    beta_agg_df <- aggregate_probes(map, beta_vals)
    ann_df <- ann[rownames(m_agg_df), , drop = FALSE]
    
    return_list[[paste0(prefix, "_", type, "_agg_m")]] <- m_agg_df
    return_list[[paste0(prefix, "_", type, "_agg_beta")]] <- beta_agg_df
    return_list[[paste0(prefix, "_", type, "_agg_ann")]] <- ann_df
  }
  
  # Extract feature-based cuts of the methylation values
  for(i in 1:length(locus_ann_list)){
    type <- names(locus_ann_list)[i]
    ann <- locus_ann_list[[i]]
    message("Cutting ", type)
    
    # Cut
    m_cut_df <- cut_probes(ann, m_vals)
    beta_cut_df <- cut_probes(ann, beta_vals)
    ann_df <- ann[rownames(m_cut_df), ]
    
    return_list[[paste0(prefix, "_", type, "_m")]] <- m_cut_df
    return_list[[paste0(prefix, "_", type, "_beta")]] <- beta_cut_df
    return_list[[paste0(prefix, "_", type, "_ann")]] <- ann_df
  }
  
  return(return_list)
}

#' Removes probes with P-value > 0.01 in any sample, overlapping SNPs
#' or annotated as cross-reactive
#' Aggregates M-values and Beta-values to regions by averaging over unfiltered loci
#' in each region. Generates annotation for resulting data
filter_and_aggregate <- function(m_set, 
                                 det_p, 
                                 cross_reactive_probes, 
                                 region_map_list, 
                                 region_ann_list,
                                 locus_ann_list,
                                 dataset_id = NA,
                                 det_p_val = 0.01,
                                 m_val_pcount = 1,
                                 prefix = "methyl"){
  ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  # Filter loci --------------------------------------------------------------------------------
  message("Filtering...")
  
  if(!all(featureNames(m_set) == rownames(det_p))) warning("orders do not match")
  keep <- rowSums(det_p < det_p_val) == ncol(m_set)
  m_set <- m_set[keep, ]
  
  # filter out probes overlapping snps
  gm_set <- mapToGenome(m_set)
  gm_set <- dropLociWithSnps(gm_set)
  
  # filter out cross-reactive probes
  keep <- !(featureNames(gm_set) %in% cross_reactive_probes$TargetID)
  gm_set <- gm_set[keep, ]
  
  # Aggregate values by region --------------------------------------------------------------------------------
  message("Aggregating and cutting M and Beta-values by region...")
  
  beta_vals <- getBeta(gm_set)
  meth_vals <- getMeth(gm_set)
  unmeth_vals <- getUnmeth(gm_set)
  m_vals <- log2((meth_vals + m_val_pcount) / (unmeth_vals + m_val_pcount))
  
  return_list <- aggregate_values(m_vals, beta_vals, region_map_list, region_ann_list, locus_ann_list, prefix)
  
  return_list[["beta_vals"]] <- beta_vals
  return_list[["m_vals"]] <- m_vals
  
  return(return_list)
}

#' Assess and summarize normalization method performance
assess_method <- function(b1_1, b1_2, b2_1, b2_2, vals, method_name){
  message("Calculating means...")
  mean1_1 <- rowMeans(vals[,b1_1], na.rm = TRUE)
  mean1_2 <- rowMeans(vals[,b1_2], na.rm = TRUE)
  mean2_1 <- rowMeans(vals[,b2_1], na.rm = TRUE)
  mean2_2 <- rowMeans(vals[,b2_2], na.rm = TRUE)
  
  message("Calculating rmse...")
  rmse1 <- sqrt(mean((mean1_2 - mean1_1)^2, na.rm = TRUE)) / iqr(mean1_2, na.rm = TRUE)
  rmse2 <- sqrt(mean((mean2_2 - mean2_1)^2, na.rm = TRUE)) / iqr(mean2_2, na.rm = TRUE)
  rmse_bw <- sqrt(mean((mean2_1 - mean1_1)^2, na.rm = TRUE)) / iqr(mean2_1, na.rm = TRUE)
  
  res <- data.frame(method = method_name,
                    stat = c("batch_1", "batch_2", "between"),
                    rmse = c(rmse1, rmse2, rmse_bw))
  
  mean_vals <- data.frame(cbind(mean1_1, mean1_2, mean2_1, mean2_2))
  
  return(list(mean_vals, res))
}