#'
#'
#'
#' @author Simon Schulte
#' Date: 2024-01-31 17:12:27.951766
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
library(disaggR)
library(arrow)
library(ggthemes)
library(testthat)
library(furrr)
library(qs)

############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())

plan(multisession, workers = 4)


path2CT <- './data/prepare_CT_UNFCCC.RData'
path2CountryMappingDESIRE <-  './data/CountryMappingDESIRE.xlsx'

path2output <- './data/'
N <- 10000

############################################################################## #
##### functions #################################################################
############################################################################## #
source('functions.R')

############################################################################## #
##### load data #############################################################
############################################################################## #
# sampling_fun <- readRDS('./intermediate_results/german_aea_sampling_fun1.RData')
# sample <- sampling_fun(N)

sample <- qread('./intermediate_results/german_aea_sample1.qs')
ct <- readRDS(path2CT)
ct <- ct[, .(CRF_code, CRF_class, EXIOBASE_code, proxy_data_source,
             EXIOBASE_products)]

# _a) data preparation =========================================================

ct[, CRF_class := normalize_UNFCCC_classfications(CRF_class)]
ct <- ct[, .(id, CRF_code, CRF_class, EXIOBASE_code, proxy_data_source, EXIOBASE_products)]

#setkey(samples, year, party, gas, category_code, classification)
#setkey(ct, CRF_code, CRF_class)



# merge sample with ct ========================================================
dt <- merge(sample, ct,
            by.x = c('category_code', 'classification'),
            by.y = c('CRF_code', 'CRF_class'),
            all.x = TRUE, sort = FALSE)
dt[, n_corres := list(sapply(EXIOBASE_code, length))]
dt[EXIOBASE_code == 'all', n_corres := 163]

# _c) checks ==================================================================

# __i. proxy data sources

# flatten list
dt$proxy_data_source %>% unique

dt[sapply(proxy_data_source, is.null),
   proxy_data_source := NA]

unlist(dt$proxy_data_source) %>% length
length(dt$proxy_data_source)


if (all(sapply(dt$proxy_data_source, function(x) length(x) <= 1))) {
  dt$proxy_data_source <- unlist(dt$proxy_data_source, recursive = FALSE)
} else {
  warning('column proxy_data_source could not be converted to character (i.e. unlisted) because at least one element has length > 1')
}
dt$proxy_data_source %>% unique

# correction 1: eurostat --> SUPPLY (only very minor emissions sources)
dt[proxy_data_source == 'eurostat',
   `:=`(proxy_data_source = 'SUPPLY',
        EXIOBASE_products = 'all')]

# correction 2: IEA: TODO insert IEA data on power-to-heat ratios
dt[proxy_data_source == 'IEA']$category_code %>% unique

# check if all information is available
dt[is.na(proxy_data_source) & n_corres > 1] # should be empty

# check for 'all' in EXIOBASE correspondences
dt[EXIOBASE_code == 'all']

# check coverage by each proxy data source:
dt$proxy_data_source %>% unique
dt[n_corres > 1 | proxy_data_source == 'PEFA',
   sum(emissions_CRF),
   by = .(gas, proxy_data_source)] %>%
  ggplot(aes(x = gas, y = V1, fill = proxy_data_source)) +
  geom_col(position = 'stack') +
  scale_fill_colorblind7() +
  facet_wrap(~gas, scales = 'free')

# check coverage ====
unfccc_samples3 <- merge(sample[, list(total_emissions = sum(emissions_CRF)),
                                by = gas],
                         dt[, list(covered_emissions=  sum(emissions_CRF)), by = gas],
                         by = 'gas')
unfccc_samples3[, percent_covered := covered_emissions / total_emissions]
unfccc_samples3 # all covered


# _d) more preparations ======================================================

dt[, key := paste(paste(EXIOBASE_code, sep = ''),
                  paste(EXIOBASE_products, sep = ''),
                  paste(proxy_data_source, sep = ''),
                  sep = '')]

# list all proxy data sources
dt$proxy_data_source %>% unique
dt$proxy_data_source %>% sapply(length) %>% unique
dt[sapply(proxy_data_source, length) == 0]
#dt[, proxy_data_source := sapply(dt$proxy_data_source, '[', 1)]

dt[is.na(proxy_data_source) & n_corres == 1
   , proxy_data_source := '1to1']
dt[proxy_data_source == 'output']
dt[proxy_data_source == 'eurostat']
dt[EXIOBASE_code == 'all']

dt[proxy_data_source == ('USE')]

dt[is.na(proxy_data_source) & n_corres != 1]

# there are still 1to1 correspondences which are declared as 'SUPPLY' --> ajust proxy data source
dt[n_corres == 1 & !(proxy_data_source %in% c('1to1', 'PEFA'))
   , proxy_data_source := '1to1']

dt[, EXIOBASE_region_code := 'DE']

############################################################################## #
##### 2. Get proxies ##########################################################
############################################################################## #

############################################################################## #
# _a) SUT ------------------------------------------------------------------
############################################################################## #

# __i. load data ---------------------------------------------------------------

use <- read_feather(file.path(path2output, 'prepare_SUT_proxies_use.feather'))
#use[, country_industry := countrycode(country_industry,'iso2c', 'iso3c')]
use <- na.omit(use)


supply <- read_feather(file.path(path2output, 'prepare_SUT_proxies_supply.feather'))
#supply <- readRDS('./temp_results/4a_supply.RData')
#supply[, country_industry := countrycode(country_industry,'iso2c', 'iso3c')]
supply <- na.omit(supply)

# load meta data
#meta <- parse_EB3_metadata('/home/simon/Documents/PhD_PROSET/data/EXIOBASE3/V3.8.2/IOT_2015_pxp')


# __ii. data checks -----------------------------------------------------------
# are all proxy details there? (should return empty data.table)

test_that('all proxy details (corresponding EB indudstries and products) are there', {
  expect_equal(0, nrow(dt[proxy_data_source == 'USE'& is.na(EXIOBASE_products)]))
  expect_equal(0, nrow(dt[proxy_data_source == 'USE'& is.na(EXIOBASE_code)]))
  expect_equal(0, nrow(dt[proxy_data_source == 'SUPPLY'& is.na(EXIOBASE_products)]))
  expect_equal(0, nrow(dt[proxy_data_source == 'SUPPLY'& is.na(EXIOBASE_code)]))
})




# all information on EB sectors correct?
#dt$EXIOBASE_code %>% unlist %>% unique %>% non_common_elements(meta$industries$CodeNr %>% unique)
#dt$EXIOBASE_products %>% unlist %>% unique %>% non_common_elements(meta$products$CodeNr %>% unique)


# __iii. retrieve information from USE table ----------------------------------
dt[proxy_data_source == 'USE'
   , proxies := list(mapply(FUN = get_SUT_shares,
                            products = EXIOBASE_products,
                            industries = EXIOBASE_code,
                            industry_countries = EXIOBASE_region_code,
                            MoreArgs = list(sut = use),
                            SIMPLIFY = FALSE)),
   by = key]

# __iv. retrieve information from SUPPLY table ---------------------------------
dt[proxy_data_source == 'SUPPLY'
   , proxies := list(mapply(FUN = get_SUT_shares,
                            products = EXIOBASE_products,
                            industries = EXIOBASE_code,
                            industry_countries = EXIOBASE_region_code,
                            MoreArgs = list(sut = supply),
                            SIMPLIFY = FALSE)),
   by = key]

# __v. check results ----------------------------------------------------------
dt
dt[sapply(proxies, function(x) if (is.data.table(x)) nrow(x) == 1 else FALSE)]

############################################################################## #
# _b) ROAD TRANSPORT ------------------------------------------------------------------
############################################################################## #
road <- read_feather(file.path(path2output, 'prepare_ROAD_TRANSPORT_proxies.feather'))
road[, proxies := lapply(proxies, as.data.table)]
road[, proxy_data_source := 'PEFA']
#road$proxies_ROAD[[1]]
setnames(road, c('region', 'proxies'), c('EXIOBASE_region_code', 'proxies_ROAD'))

dt <- merge(dt, road,
            by = c('EXIOBASE_region_code', 'proxy_data_source',
                   'gas'),
            sort = FALSE,
            all.x = TRUE)

dt[sapply(proxies, is.data.table) & sapply(proxies_ROAD, is.data.table)]
dt[sapply(proxies_ROAD, is.data.table), proxies := proxies_ROAD]
dt[, proxies_ROAD := NULL]




############################################################################## #
# _c) IEA (power-to-heat ratio) ------------------------------------------------------------------
############################################################################## #
# atm: use 0.5 (average of default power-to-heat values of table 2 of below report)
# see: https://ec.europa.eu/eurostat/documents/38154/42195/Final_CHP_reporting_instructions_reference_year_2016_onwards_30052017.pdf/f114b673-aef3-499b-bf38-f58998b40fe6
dt[proxy_data_source == 'USE']$proxies[[1]]

dt[proxy_data_source == 'IEA',
   proxies := list(list(data.table(
     industry_code = unlist(EXIOBASE_code, recursive = FALSE),
     share = 1 / length(unlist(EXIOBASE_code))
   ))),
   by = .(key, gas)]


############################################################################## #
##### Tests #################################################################
############################################################################## #

shares_not_sum_to_1 <-   dt[
  !is.na(proxy_data_source)
  & !(proxy_data_source %in% c('1to1', 'PEFA'))
  & !sapply(proxies, function(x) isTRUE(all.equal(sum(x$share), 1)))
]

test_that('all proxy shares sum to 1 (SUPPLY, USE, IEA)', {
  expect_equal(nrow(shares_not_sum_to_1), 0)
})

missing_proxies <- dt[!(proxy_data_source %in% c('1to1'))
                      & !sapply(proxies, is.data.table)
                      , .(category_code, classification)] %>%
  unique


test_that('no proxies are missing', {
  expect_equal(nrow(missing_proxies), 0)
})


############################################################################## #
##### prepare data for sampling #################################################################
############################################################################## #

sample_cols <- grep('sample_case_', names(dt), value = TRUE)
dt <- dt[, -c('key', 'n_corres', 'EXIOBASE_products', 'id')]

#dt <- dt[, .(ID, gas, category_code, classification, EXIOBASE_region_code,
#             emissions_CRF, sample_1, sample_2, EXIOBASE_code, proxies, proxy_data_source)]

# extract alpha
dt[, alpha2 := lapply(proxies, function(x) {
  setNames(x$share, x$industry_code)
})]
dt[proxy_data_source == '1to1',
   alpha2 := pmap(list(x = EXIOBASE_code), function(x) setNames(1, x))]

# create uninformative alpha
dt[, alpha1 := pmap(list(x = alpha2), function(x) setNames(rep(1/length(x),
                                                               length(x)),
                                                           names(x)))]

# make up beta
dt[, beta1 := pmap(list(x = alpha2), function(x) setNames(rep(NA, length(x)), names(x)))]

dt[, beta2 := beta1]
dt[proxy_data_source == 'USE',
   beta2 := pmap(list(x = alpha2), function(x) setNames(0.1 * x, names(x)))]
dt[proxy_data_source == 'SUPPLY',
   beta2 := pmap(list(x = alpha2), function(x) setNames(0.2 * x, names(x)))]
dt[proxy_data_source == 'PEFA',
   beta2 := pmap(list(x = alpha2), function(x) setNames(0.3 * x, names(x)))]


dt[proxy_data_source == 'USE',
   beta3 := pmap(list(x = alpha2), function(x) setNames(0.005 * x, names(x)))]
dt[proxy_data_source == 'SUPPLY',
   beta3 := pmap(list(x = alpha2), function(x) setNames(0.01 * x, names(x)))]
dt[proxy_data_source == 'PEFA',
   beta3 := pmap(list(x = alpha2), function(x) setNames(0.015 * x, names(x)))]


dt[, .(beta1, beta2, beta3)] %>%
  .[, lapply(.SD, unlist)] %>%
  melt %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = variable),
                 alpha = 0.3, position = 'identity') +
  # geom_boxplot(aes(x = variable, y = value)) +
  scale_x_log10() +
  scale_fill_colorblind()

dt[, .(alpha1, alpha2)] %>%
  .[, lapply(.SD, unlist)] %>%
  melt %>%
  ggplot() +
  geom_histogram(aes(x = value, fill = variable),
                 alpha = 0.3, position = 'identity') +
  # geom_boxplot(aes(x = variable, y = value)) +
  scale_x_log10() +
  scale_fill_colorblind()




#
# #
# test <- dt[proxy_data_source == 'PEFA' & gas == 'CO2'][1]
# qsave(test, 'test_data_extreme.qs')

############################################################################## #
##### Sample #############################################################
############################################################################## #

sample_cols <- grep('sample_case_', names(dt), value = TRUE)
dt2sample <- dt[, c('ID','gas', 'alpha1', 'alpha2',
       'beta1', 'beta2', 'beta3', 'EXIOBASE_code', sample_cols), with = FALSE]
 # dt[, .(ID, gas, sample_1, sample_2, alpha1, alpha2,
 #                    beta1, beta2, EXIOBASE_code)]

# dt2sample$sample_1 <- as.list(dt2sample$sample_1)
# dt2sample$sample_2 <- as.list(dt2sample$sample_2)


# select only co2 emissions (change for later)
dt2sample <- dt2sample[gas == 'CO2']

# make up different combinations
cases <- as.data.table(tribble(
  ~case, ~x, ~alpha, ~beta,
  1, 'sample_case_1', 'alpha1', 'beta1',
  2, 'sample_case_2', 'alpha1', 'beta1',
  3, 'sample_case_2', 'alpha2', 'beta1',
  4, 'sample_case_4', 'alpha2', 'beta2',
  5, 'sample_case_4', 'alpha2', 'beta3'
))

# sample
cases[, sample := furrr::future_pmap(list(x = x, alpha = alpha, beta = beta),
                       function(x, alpha, beta) {
                         sample_system(N, data = dt2sample, x = x, shares = alpha,
                                       sds = beta, shares_lb = 1E-2)
                       }, .progress = TRUE, .options = (furrr_options(seed = TRUE)))]

# convert matrices to data.table
cases[, sample := pmap(list(x = sample), function(x) {
  x <- convert_sample_to_dt(x)
  setkey(x, y)
  return(x)
})]

############################################################################## #
##### save results #############################################################
############################################################################## #
qsave(cases, './intermediate_results/german_aea_sample2.qs')





# THE END ---------------------------------------------------------------------
