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

# library(WGCNA)
library(data.table)
library(tidyverse)
library(units)
library(ggforce)
library(qs)
library(mRio)
library(Matrix)
library(furrr)
library(ggthemes)

############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())

n_cores <- 6

############################################################################## #
##### functions #################################################################
############################################################################## #
source('functions.R')

#plan(multisession, workers = 4)


############################################################################## #
##### load data #############################################################
############################################################################## #
fp_multiplier <- qread('./intermediate_results/german_aea_fp_multiplier.qs')
aea_raw <- qread('./intermediate_results/german_aea_sample2.qs')



# calc correlations
fp_multiplier2 <- fp_multiplier[, list(data = list(data.table(i, j, value))),
                                by = case]
fp_multiplier2[, data := pmap(list(x = data), function(x) dcast(x, i ~j,
                                                                value.var = 'value'))]
fp_multiplier2[, cor := pmap(list(x = data), function(x) cordt(x[, -'i'],
                                                               nThreads = n_cores,
                                                               large = FALSE))]

fpm_cor_sample <- fp_multiplier2[, rbindlist(cor), by = case]



# accounts
aea <- aea_raw[, rbindlist(sample), by = case]
N <- aea$sample[[1]] %>% length
aea_summary <- calculate_summary_statistics(aea)
aea_sample <- aea[, list(value = unlist(sample, FALSE)), by = .(case, y)]

aea_sample[, run_id := rep(1:N, nrow(aea))]



aea_cor <- aea_sample[, list(data = list(data.table(run_id, y, value))),
                      by = case]

aea_cor[, data := pmap(list(x = data), function(x) dcast(x, run_id ~y,
                                                         value.var = 'value'))]
aea_cor[, cor := pmap(list(x = data), function(x) cordt(x[, -'run_id'],
                                                        nThreads = n_cores,
                                                        use_names = TRUE,
                                                        large = FALSE))]
aea_cor <- aea_cor[, rbindlist(cor), by = case]





############################################################################## #
##### save results #############################################################
############################################################################## #

qsave(fpm_cor_sample, './intermediate_results/german_aea_fp_multiplier_cor.qs')
qsave(aea_cor, './intermediate_results/german_aea_cor.qs')


