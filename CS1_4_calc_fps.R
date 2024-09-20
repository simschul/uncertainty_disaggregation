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
library(MethylCapSig)
library(proxyC)
############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())

path2output <- "./data"
path2exiobaseIOT <- "/home/simon/Documents/PhD_PROSET/data/EXIOBASE3/V3.8.2/IOT_2015_ixi"
############################################################################## #
##### functions #################################################################
############################################################################## #
source('functions.R')

#plan(multisession, workers = 4)
############################################################################## #
##### load data #############################################################
############################################################################## #
# load samples
dt_raw <- qread('./intermediate_results/german_aea_Fmat.qs')

dt_raw[, value := lapply(value, matsbyname::Matrix)]


indices_F <- qread('./intermediate_results/german_aea_indices_F.qs')
dt_mats <- dt_raw
############################################################################## #
##### load EXIOBASE #############################################################
############################################################################## #
xvec <- readRDS(file.path(path2output, 'prepare_EXIOBASE_x.RData'))
Ymat <- readRDS(file.path(path2output, 'prepare_EXIOBASE_Y.RData'))
Lmat <- readRDS(file.path(path2output, 'prepare_EXIOBASE_L.RData'))

# prepare L
colnames(Lmat) <- as.character(1:ncol(Lmat))
rownames(Lmat) <- as.character(1:nrow(Lmat))
Lmat <- Lmat[indices, indices]
Lmat <- matsbyname::Matrix(Lmat)
dim(Lmat)


# prepare x
names(xvec) <- as.character(1:length(xvec))
xvec <- xvec[indices]
xvec1 <- 1/xvec
xvec1[!is.finite(xvec1)] <- 0
#xvec_sparse <- sparseVector(xvec1, 1:length(xvec1), length = length(xvec1))
xhat <- Diagonal(n = length(xvec),x = xvec1)
colnames(xhat) <- indices
rownames(xhat) <- indices
xhat <- matsbyname::Matrix(xhat)

# prepare Y
rownames(Ymat) <- as.character(1:nrow(Ymat))
Ymat2 <- matsbyname::Matrix(Ymat)
Ymat2 <- Ymat2[indices, 'DE', drop = FALSE]
dim(Ymat2)


# calculate S
dt_mats[, Smat := pmap(list(x = value), function(x) x %*% xhat)]

# calculate footprints
dt_mats[, fp_multiplier := pmap(list(x = Smat), function(x){
  res <- x %*% Lmat
  #res[, -ids_DE] <- 0
  return(res)
} ,.progress = TRUE)] # 20 sec for 100 runs x 4 cases
dt_mats[, fp_national := pmap(list(x = fp_multiplier), function(x) {
  x %*% Ymat2
} ,.progress = TRUE)]


# convert to data.tab.e
dt_mats[, fp_multiplier_dt := pmap(list(x = fp_multiplier), function(x) {
  convert_matsbyname_to_dt(x)
})]
dt_mats[, fp_national_dt := pmap(list(x = fp_national), function(x) {
  convert_matsbyname_to_dt(x)
})]


# extract data
fp_multiplier <- dt_mats[variable == 'Fmat', rbindlist(fp_multiplier_dt), by = .(case)]
fp_national <- dt_mats[variable == 'Fmat', rbindlist(fp_national_dt), by = .(case)]



fp_multiplier_resampled <- dt_mats[,rbindlist(fp_multiplier_dt),
                                   by = .(case, variable)]
fp_national_resampled <- dt_mats[,rbindlist(fp_national_dt),
                                   by = .(case, variable)]



############################################################################## #
##### save results #############################################################
############################################################################## #
qsave(fp_multiplier, './intermediate_results/german_aea_fp_multiplier.qs')
qsave(fp_national, './intermediate_results/german_aea_fp_national.qs')

qsave(fp_multiplier_resampled, './intermediate_results/german_aea_fp_multiplier_resampled.qs')
qsave(fp_national_resampled, './intermediate_results/german_aea_fp_national_resampled.qs')
#qsave(fp_multiplier2[, .(case, cor)], './intermediate_results/german_aea_fp_multiplier_cor.qs')

# THE END ---------------------------------------------------------------------
