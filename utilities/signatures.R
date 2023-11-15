#' Assign types based on the output of local.gsl.permutation
#' function from hexplot analysis. Assign none if no type
#' has probability > min_prob. Otherwise, assign to the type 
#' or combination of types with probability > min_prob
prob2type <- function(prob_data, min_prob = 0.95){
  type_ids <- apply(prob_data, 1, function(x){
    flag_vec <- x > min_prob
    if(all(flag_vec)){
      "all"
    } else {
      paste0(colnames(prob_data)[flag_vec], collapse = "_and_")
    }
  })
  return(type_ids)
}

#' perform permutation test on expression matrix and genesets structures
#' exp_matrix: expression matrix with samples in columns
#' gsl: list of signatures with format of -- gsl = list(sample1=list(up=(), dn=()), sample2=list(up=(), dn=()), sample3=list(up=(), dn=()))
score_gsl <- function(exp_matrix, gsl, N=1000, Nthreads=10){
  # load library
  require(foreach)
  require(doMC)
  registerDoMC(Nthreads)
  # calculate stats from test geneset
  test.stats = sapply(gsl, function(x){
    colMeans(exp_matrix[x$up,])-colMeans(exp_matrix[x$dn,])
  })
  # init progress
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  # get permutation stats for N times
  rand.test = foreach(i = 1:N) %dopar% {
    # make geneset based on the structure of gsl
    allID = unique(unlist(gsl))
    rand.N = length(allID)
    set.seed(i*12345)
    rand.gene = sample(rownames(exp_matrix), rand.N)
    rand.gsl = lapply(gsl, function(x){
      lapply(x, function(y){
        rand.gene[allID %in% y]
      })
    })
    # calculate stats
    rand.stats = sapply(rand.gsl, function(x){
      colMeans(exp_matrix[x$up,])-colMeans(exp_matrix[x$dn,])
    })
    # update progress
    setTxtProgressBar(pb, i)
    # compare to test stats
    test.stats-rand.stats
  }
  # permutation pvalue
  perm.pval = Reduce('+', rand.test)/N
  # output
  data.frame(perm.pval)
}

#' Assign types based on the output of easy_signature_score
#' Assert there are no ties for top score, and then label
#' as the top-scoring type, with labels inferred from 
#' column names of input data
score2type <- function(score_data){
  n_top <- apply(score_data, 1, function(x){
    length(which(x == max(x)))
  })
  n_ties <- length(which(n_top > 1))
  stopifnot(n_ties == 0)
  assigned_types <- apply(score_data, 1, function(x){
    colnames(score_data)[which.max(x)]
  })
}

#' Calculate simple signature scores for individual samples for the given
#' gsl
easy_signature_score <- function(exp_data, gsl){
  sapply(gsl, function(x){
    colMeans(exp_data[x$up, ]) - colMeans(exp_data[x$dn, ])
  })
}

#' Average expression per sample
easy_gene_set_score <- function(exp_data, gsl){
  sapply(gsl, function(x){
    colMeans(exp_data[x$gene,])
  })
}
