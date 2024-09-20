#'
#'
#'
#' @author Simon Schulte
#' Date: 2024-02-07 09:40:31.284218
#'
#' Content:
#'


############################################################################## #
##### load packages ############################################################
############################################################################## #

library(data.table)
library(tidyverse)
library(units)
library(ggforce)
library(qs)
library(mRio)
library(Matrix)
library(MatrixExtra)
library(sparseMatrixStats)
library(furrr)
library(ggthemes)
library(WGCNA)
# library(MethylCapSig)
library(truncnorm)


############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())

path2exiobaseIOT <- "/home/simon/Documents/PhD_PROSET/data/EXIOBASE3/V3.8.2/IOT_2015_ixi"
############################################################################## #
##### functions #################################################################
############################################################################## #
source('functions.R')
source('function_lcmix.R')
#plan(multisession, workers = 4)


############################################################################## #
##### load data #############################################################
############################################################################## #

# load samples
dt_raw <- qread('./intermediate_results/german_aea_sample2.qs')


N <- length(dt_raw[['sample']][[1]]$sample[[1]])

# prepare sample data.tables to merge
dt_raw[, sample := pmap(list(x = sample), function(x) {
  x[, region := 'DE']
  x[, row := 1]
  setnames(x, 'y', 'CodeNr')
  return(x)
})]

# load meta data
meta <- parse_EB3_metadata(path2exiobaseIOT)
indices_S <- read_EB3_S_meta(file.path(path2exiobaseIOT, 'satellite', 'F.txt'))
indices_S$colnames[, col := 1:.N]
indices_S$rownames[, row := 1:.N]
indices_S$colnames <- merge(indices_S$colnames, meta$industries, by.x = 'sector',
                            by.y = 'Name',
                            sort = FALSE)
ids_DE <- indices_S$colnames[region == 'DE']$col

# merge both
dt_raw[, sample := pmap(list(x = sample), function(x) {
  merge(x, indices_S$colnames[, .(region, col, CodeNr)],
        by = c('region', 'CodeNr'), all.x = TRUE)
})]

#dt_raw[, sample := pmap(list(x= sample), function(x) x[, row := 1])]

# convert to matrix
dt_raw[, Fmat := pmap(list(x = sample), function(x) {
  x <- x[, unlist_sample(sample), by = col]
  x <- na.omit(x)
  setorder(x, col, id_run)
  mat <- data.table::dcast(x, id_run ~ col, value.var = "value", fill = 0)
  rownames <- mat[['id_run']]
  mat <- as(mat[, -'id_run'], 'sparseMatrix')
  rownames(mat) <- rownames
  return(mat)
  # sparseMatrix(x[['id_run']], x[['col']], x = x[['value']],
  #              dims = c(N, 7987))
})]


# dt_raw[, Fmat := pmap(list(x = sample), function(x) {
#   as_sparse_matrix_list(na.omit(x), nrow = 1, ncol = 7987)
# }, .progress = TRUE)]



#dt_mats <- dt_raw[, list(Fmat = list_flatten(Fmat)), by = case]
#dt_mats[, run_id := rep(1:N, length(unique(dt_raw$case)))]

# dt_mats <- dt_mats[run_id < 1001] # TODO: include all samples, delete this line


# resample F matrices (using both uni- and multivariate lognorm)
dt_raw[, means := pmap(list(x = Fmat), function(x) {
  out <- Matrix::colMeans(x, na.rm = TRUE)
  #as(out, 'sparseVector')
  out
})]
dt_raw[, means_log := pmap(list(x = Fmat), function(x) {
  x@x <- log(x@x)
  out <- Matrix::colMeans(x, na.rm = TRUE)
  # as(out, 'sparseVector')
  out
})]
dt_raw[, vars := pmap(list(x = Fmat), function(x) {
  sparseMatrixStats::colVars(x, na.rm = TRUE)
  #x <- apply(x, 2, var)
  #x <- colVars(x)
  #ifelse(!is.finite(x) , 0, x)
})]
dt_raw[, sds := pmap(list(x = Fmat), function(x) {
  sparseMatrixStats::colSds(x, na.rm = TRUE)
  # x <- proxyC::colSds(x)
  # ifelse(!is.finite(x) , 0, x)
})]
dt_raw[, sds_log := pmap(list(x = Fmat), function(x) {
  x@x <- log(x@x)
  sparseMatrixStats::colSds(x, na.rm = TRUE)
  #x <- proxyC::colSds(x)
  #x <- apply(x, 2, sd)
  #ifelse(!is.finite(x) , 0, x)
})]
dt_raw[, cor_log := pmap(list(x = Fmat), function(x) {
  x@x <- log(x@x)
  #x <- WGCNA::cor(x)
  x <- sparse_cor(x)
  x <- Matrix(x, sparse = TRUE)
  return(x)
  # x <- proxyC::simil(x, method = 'correlation')
})]

dt_raw[, cor := pmap(list(x = Fmat), function(x) {
  #x <- WGCNA::cor(x)
  x <- sparse_cor(x)
  #x <- Matrix(x, sparse = FALSE)
  return(x)
  # x <- proxyC::simil(x, method = 'correlation')
})]


dt_raw[, cov := pmap(list(x = Fmat), function(x) {
  x <- sparse_cov(x)
  return(x)
})]





# resample
dt_raw[, Fmat_resampled_mv := pmap(list(n = N, Mu = means, Sigma = vars,
                                        R = cor_log), mvlognormal)]
dt_raw[, Fmat_resampled_mv := pmap(list(x = Fmat_resampled_mv, y = Fmat), function(x, y) {
  x <- Matrix(x, sparse = TRUE)
  colnames(x) <- colnames(y)
  rownames(x) <- rownames(y)
  return(x)
})]

# dt_raw[, Fmat_resampled_uv := pmap(list(n = N, meanlog = means_log,
#                                         sdlog = sds_log), rlnorm2)]
dt_raw[, Fmat_resampled_uv := pmap(list(n = N, mean = means,
                                        sd = sds), rlnorm2_vec)]

dt_raw[, Fmat_resampled_uv := pmap(list(x = Fmat_resampled_uv, y = Fmat), function(x, y) {
  x <- Matrix(x, sparse = TRUE)
  colnames(x) <- colnames(y)
  rownames(x) <- rownames(y)
  return(x)
})]


# gamma distribution resample
dt_raw[, Fmat_resampled_mv_gamma := pmap(list(n = N, mean = means, sd = sds,
                                        corr = cor), rmvgamma2)]
dt_raw[, Fmat_resampled_mv_gamma := pmap(list(x = Fmat_resampled_mv_gamma, y = Fmat), function(x, y) {
  x <- Matrix(x, sparse = TRUE)
  colnames(x) <- colnames(y)
  rownames(x) <- rownames(y)
  return(x)
})]

# dt_raw[, Fmat_resampled_uv := pmap(list(n = N, meanlog = means_log,
#                                         sdlog = sds_log), rlnorm2)]
dt_raw[, Fmat_resampled_uv_gamma := pmap(list(n = N, mean = means,
                                        sd = sds), rgamma2_vec)]

dt_raw[, Fmat_resampled_uv_gamma := pmap(list(x = Fmat_resampled_uv_gamma, y = Fmat), function(x, y) {
  x <- Matrix(x, sparse = TRUE)
  colnames(x) <- colnames(y)
  rownames(x) <- rownames(y)
  return(x)
})]


# trunc normal distribution resample
k <- length(dt_raw$means[[1]])
lower <- rep(0, k)

system.time({
  dt_raw[, Fmat_resampled_mv_truncnorm := pmap(list(n = N, mean = means,
                                                    sigma = cov),
                                               function(n, mean, sigma, lower = rep(0, k)) {
                                                 tmvtnorm::rtmvnorm(n = n, mean = mean, sigma= sigma,
                                                                    lower = lower,
                                                                    algorithm = 'gibbs')
                                               })]

})



dt_raw[, Fmat_resampled_mv_truncnorm := pmap(list(x = Fmat_resampled_mv_truncnorm, y = Fmat),
                                             function(x, y) {
  x <- Matrix(x, sparse = TRUE)
  colnames(x) <- colnames(y)
  rownames(x) <- rownames(y)
  return(x)
})]

dt_raw[, Fmat_resampled_uv_truncnorm := pmap(list(n = N, mean = means,
                                                  a = 0,
                                              sd = sds), rtruncnorm_vec)]

dt_raw[, Fmat_resampled_uv_truncnorm := pmap(list(x = Fmat_resampled_uv_truncnorm, y = Fmat),
                                         function(x, y) {
  x <- Matrix(x, sparse = TRUE)
  colnames(x) <- colnames(y)
  rownames(x) <- rownames(y)
  return(x)
})]

dt_raw <- melt(dt_raw[, .(case, Fmat, Fmat_resampled_mv, Fmat_resampled_uv,
                          Fmat_resampled_mv_gamma, Fmat_resampled_uv_gamma,
                          Fmat_resampled_mv_truncnorm, Fmat_resampled_uv_truncnorm)],
     id.vars = c('case'))


############################################################################## #
##### save results #############################################################
############################################################################## #
qsave(dt_raw, './intermediate_results/german_aea_Fmat.qs')

# qsave(dt_raw[, .(case, Fmat, Fmat_resampled_mv, Fmat_resampled_uv, Fmat_resampled_mv_gamma,
#                  Fmat_resampled_uv_gamma)],
#       './intermediate_results/german_aea_Fmat.qs')
qsave(indices_S$colnames, './intermediate_results/german_aea_indices_F.qs')


# THE END ---------------------------------------------------------------------
