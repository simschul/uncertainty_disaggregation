#'
#'
#'
#' @author Simon Schulte
#' Date: 2024-01-04 13:24:45.238351
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
library(data.tree)
library(collapsibleTree)
library(testthat)
library(logr)
library(disaggR)
library(arrow)
library(furrr)
library(qs)


############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())

#plan(multisession, workers = 6)

path2paper2results <- './intermediate_results/'

COUNTRY <- 'DEU'
YEAR <- '2015'
GAS <- 'CO2'
N <- 10000
############################################################################## #
##### functions #################################################################
############################################################################## #
source(file.path('functions.R'))
source(file.path('src', 'functions_dirichlet.R'))
source(file.path('src', 'functions_trees.R'))


############################################################################## #
##### load data #############################################################
############################################################################## #

# load data
dt <- readRDS(file.path(path2paper2results, 'merge_CRF_NIR_data.RData'))
NATIONAL_TOTALS <- readRDS(file.path(path2paper2results, 'prepare_UNFCCC_CRF_NationalTotals.RData'))

### workaround: handle missing uncertainty estimates
### calculate average CV for each CRF category
dt_cvmean <- dt[, list(cv_DEFAULT = mean(cv_NIR, na.rm = TRUE)),
                by = .(gas, category_code)] %>%
  na.omit

dt <- merge(dt, dt_cvmean, by = c('gas', 'category_code'),
            all.x = TRUE, sort = FALSE)
dt[is.na(cv_DEFAULT), cv_DEFAULT := median(dt_cvmean$cv_DEFAULT)]
dt[, sd_DEFAULT := cv_DEFAULT * emissions_CRF]

# select
dt <- dt[party == COUNTRY & year == YEAR & gas == GAS]

#### begin
dt[, ID := 1:.N]
dt[, sd_CRF := emissions_CRF * cv_NIR]

# subset data.table: only nodes that have unvertainty information and are no leaves
dt_full_nodes <- dt[!is.na(cv_NIR), #& is.na(is_leaf_leaf),
                    .(ID, gas, id, classification_hierarchy, emissions_CRF,
                      sd_CRF, sd_DEFAULT, cv_DEFAULT)]

NATIONAL_TOTALS[party == COUNTRY & year == YEAR]
dt_full_nodes[, sum(emissions_CRF), by = gas]


for (i in 1:nrow(dt_full_nodes)) {
  irow <- dt_full_nodes[i, .(ID, gas, id,
                             classification_hierarchy,
                             emissions_CRF, cv_DEFAULT)] %>%
    as.list

  dt_subcats <- dt[id %in% find_subcategories(id, irow$id) & gas == irow$gas &
                     is_leaf_leaf == TRUE]

  if (nrow(dt_subcats) > 0) {
    # node has subcategories
    test_that('emissions of subcategory are coherent', {
      expect_equal(dt_subcats[, sum(emissions_CRF)], irow$emissions_CRF)
    })

    # calculate proxy data
    dt_subcats[, proxy := emissions_CRF / sum(emissions_CRF)]
    dt_subcats[, beta := proxy * cv_DEFAULT]
    dt_subcats[, .(ID, proxy, beta)]

    # add to data.table
    dt_full_nodes[ID == irow$ID,
                  `:=`(alpha = list(dt_subcats$proxy),
                       beta = list(dt_subcats$beta),
                       y = list(dt_subcats$ID)
                  )]
  } else {
    # node has no subcategories
    dt_full_nodes[ID == irow$ID,
                  `:=`(alpha = list(1),
                       beta = list(NA),
                       y = list(irow$ID)
                  )]
  }


  cat(i, '')
}

# check that all alphas sum to one
test_that('all alphas sum to one', {
  expect_equal(0, nrow(dt_full_nodes[sapply(alpha, function(x) !isTRUE(all.equal(sum(x), 1))), ]))
})

# are there any leaves not invluded? (neither as source, nor target)

IDs_incl <- unique(c(
  dt_full_nodes$ID, # source
  unlist(dt_full_nodes$y), # target
  unique(unlist(lapply(dt_full_nodes$id, find_parent_categories,
                       all_codes = dt$id))) # parent categories
))
dt_missing <- dt[is_leaf_leaf == TRUE & !(ID %in% IDs_incl)]
dt_missing[, alpha := list(list(1))]
dt_missing[, beta := list(list(NA))]
dt_missing[, y := lapply(ID, function(x) x)]


# prepare data to sample
#dt2sample <- dt_full_nodes[, .(ID, emissions_CRF, sd_CRF, alpha, y)]
dt2sample <- rbindlist(list('truncnorm' = dt_full_nodes[, .(ID, gas, emissions_CRF,
                                                            sd_CRF, sd_DEFAULT,
                                                            cv_DEFAULT,
                                                            alpha, beta, y)],
                            'exp' = dt_missing[, .(ID,gas, emissions_CRF, sd_CRF,
                                                   sd_DEFAULT, cv_DEFAULT, alpha, beta,
                                                   y)]),
                       idcol = 'dist')

# test if all emissions are covered now
dt2sample[, sum(emissions_CRF), by = gas]
NATIONAL_TOTALS[party == COUNTRY & year == YEAR]


dt2sample[, x0a := pmap(list(mean = emissions_CRF, sd = sd_CRF),
                        function(mean, sd) list(mean = mean, sd = sd, min = 0))]
dt2sample[, x0b := pmap(list(mean = emissions_CRF, sd = sd_DEFAULT),
                        function(mean, sd) list(mean = mean, sd = sd, min = 0))]


#TODO: add case where sd is taken to Default value (instaed of NA) to sample
# from truncnorm instead of exp

dt2sample[, beta1 := lapply(alpha, function(x) rep(NA, length(x)))]
setnames(dt2sample, 'beta', 'beta2')


dt2sample[, alpha := Map(f = function(x, y) {
  names(x) <- as.character(y)
  return(x)
} ,
x = alpha,
y = y)]

dt2sample[, beta1 := Map(f = function(x, y) {
  names(x) <- as.character(y)
  return(x)
} ,
x = beta1,
y = y)]

dt2sample[, beta2 := Map(f = function(x, y) {
  names(x) <- as.character(y)
  return(x)
} ,
x = beta2,
y = y)]


# sample
# make up different combinations
sample_dt <- as.data.table(tribble(
  ~case, ~x, ~alpha, ~beta,
  1, 'x0a', 'alpha', 'beta1',
  2, 'x0b', 'alpha', 'beta1',
  4, 'x0b', 'alpha', 'beta2'
))

# sample
sample_dt[, sample := pmap(list(x = x, alpha = alpha, beta = beta),
                       function(x, alpha, beta) {
                         sample_system(N, data = dt2sample, x = x, shares = alpha,
                                       sds = beta, shares_lb = 1E-2)
                       }, .progress = TRUE)]

# convert matrices to data.table
sample_dt[, sample := pmap(list(x = sample), function(x) {
  x <- convert_sample_to_dt(x)
  setkey(x, y)
  return(x)
})]

sample_dt2 <- sample_dt[, rbindlist(sample), by = .(case)]
sample_dt2[, case := paste0('sample_case_', case)]
sample_dt2 <- dcast(sample_dt2, y ~ case, value.var = 'sample')


# # case 1: exp when no uncertainty data is avail
# sample1 <- sample_system(N, data = dt2sample, x = 'x0a', alpha = 'alpha',
#               beta = 'beta', y = 'y')
#
# sample_dt1 <- convert_sample_to_dt(sample1)
# sampling_fun1 <- generate_sampling_fun2(data = dt2sample, x = 'x0a',
#                                        alpha = 'alpha',
#                                        beta = 'beta', y = 'y')
#
# # case 2: take default uncertainty data and sample from truncnorm
# sample2 <- sample_system(N, data = dt2sample, x = 'x0b', alpha = 'alpha',
#                          beta = 'beta', y = 'y')
#
# sample_dt2 <- convert_sample_to_dt(sample2)
# sampling_fun2 <- generate_sampling_fun2(data = dt2sample, x = 'x0b',
#                                         alpha = 'alpha',
#                                         beta = 'beta', y = 'y')
#
#
# # merge both
# sample_dt <- merge(sample_dt1, sample_dt2, by = 'y', all = TRUE, sort = FALSE,
#                    suffixes = c('_1', '_2'))
#
# add meta data
setnames(sample_dt2, 'y', 'ID')
sample_dt2[, ID := as.numeric(ID)]


sample_dt_full <- merge(sample_dt2, dt[, .(ID, gas, category_code, classification,
                                      emissions_CRF)],
      all.x = TRUE, sort = FALSE)
setkey(sample_dt_full, ID)


############################################################################## #
##### save results #############################################################
############################################################################## #
#saveRDS(sampling_fun, './intermediate_results/german_aea_sampling_fun1.RData')
qsave(sample_dt_full, './intermediate_results/german_aea_sample1.qs')


# THE END ---------------------------------------------------------------------
















