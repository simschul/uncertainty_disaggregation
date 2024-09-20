#'
#'
#'
#' @author Simon Schulte
#' Date: 2023-11-17 08:45:53.597817
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
library(gtools)
library(truncnorm)
library(mvtnorm)
# library(tmvtnorm)
library(faux)
# library(compositions)
library(disaggR)
library(ggthemes)
library(TruncatedNormal)
# library(moments)
# library(mvLognCorrEst)
library(MethylCapSig)
library(matrixStats)




############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())


theme_border <- theme_gray() +
  theme(plot.background = element_rect(fill = NA, colour = 'grey', linewidth = 1))

theme_simschul <- function(){
  #font <- "Georgia"   #assign font family up front

  theme_bw() %+replace%    #replace elements we want to change

    theme(
      strip.background  = element_rect(fill = 'grey90', colour = NA)

    )
}
theme_set(theme_simschul())
#RhpcBLASctl::blas_set_num_threads(6)
############################################################################## #
##### functions #################################################################
############################################################################## #
source('functions_dirichlet.R')
source('functions.R')
source('function_lcmix.R')


summary_statistics <- function(x) {
  temp <- list(
    min = min(x),
    max = max(x),
    mean = mean(x),
    sd = sd(x),
    Q.25 = quantile(x, probs = 0.25),
    Q.75 = quantile(x, probs = 0.75),
    Q.025 = quantile(x, probs = 0.025),
    Q.975 = quantile(x, probs = 0.075)
  )
  data.table(measure = names(temp),
             value = unlist(temp))
}

############################################################################## #
##### load data #############################################################
############################################################################## #


# 1. Set parameters ===============================================================
N <- 1E5
param_dt <- as.data.table(tribble(
  ~id_case, ~mean_0, ~sd_0, ~shares, ~sds,
  1, 100, 5, c(0.1, 0.3, 0.6), NULL,
  2, 100, 80, c(1/3, 1/3, 1/3), NULL,
  #2, 100, 10, c(0.1, 0.3, 0.6), c(0.1, 0.3, 0.6)*0.3,
  3, 100, 80, c(0.1, 0.3, 0.6), c(0.1, 0.3, 0.6)*0.05
))


# 2. Sample ===============================================================

param_dt[, sample := pmap(list(mean_0 = mean_0, sd_0 = sd_0,
                               shares = shares, sds = sds, min = 0, n= N),
                          rdisagg)]

param_dt$sample[[1]] %>% rowSums()
param_dt[, cor := pmap(list(x = sample), cordt)]

param_dt[, .(id_case, mean_0, sd_0, shares, sds)]
param_dt$cor


# 3. Resample ===============================================================
dists_dt <- as.data.table(tribble(
  ~id_dist, ~ dist, ~ type, ~name,
  1, 'lnorm', 'univariate2', 'uv_lnorm',
  2, 'lnorm', 'mvlognormal', 'mv_lnorm',
  3, 'gamma', NULL, 'uv_gamma',
  4, 'mvgamma', NULL, 'mv_gamma'
))


dt <- expand.grid(id_case = param_dt$id_case, id_dist = dists_dt$id_dist) %>%
  as.data.table %>%
  merge(param_dt, by = 'id_case') %>%
  merge(dists_dt, by = 'id_dist')

dt[, sample_new := pmap(list(dist = dist, type = type, sample = sample),
                        reconstruct_disagg_sample)]

dt[, sample_agg := lapply(sample, rowSums)]
dt[, sample_agg_new := lapply(sample_new, rowSums)]

# 4. Unnest =================================
dt2 <- dt[, .(id_dist, id_case, name, sample_agg, sample_agg_new)] %>%
  .[, list(sample = unlist(sample_agg, recursive = FALSE),
           sample_new = unlist(sample_agg_new, recursive = FALSE)),
    by = .(id_dist, id_case, name)]

dt2 <- melt(dt2, id.vars = c('id_dist', 'id_case', 'name'))

# _a) Gamma dist ======================================
dt_gamma <-  dt2[name %in% c('uv_gamma', 'mv_gamma')]

dt_gamma <- dt_gamma[!(name == 'uv_gamma' & variable == 'sample')]
dt_gamma[variable == 'sample', type := 'original']
dt_gamma[variable == 'sample_new' & name == 'mv_gamma', type := 'multivariate']
dt_gamma[variable == 'sample_new' & name == 'uv_gamma', type := 'univariate']
dt_gamma[type == 'original', alpha := 0.6]
dt_gamma[type != 'original', alpha := 0.3]
dt_gamma[, type := factor(type, levels = c('original', 'multivariate', 'univariate'))]



dt_gamma[id_case == 1, case := 'negative correlations']
dt_gamma[id_case == 2, case := 'no correlations']
dt_gamma[id_case == 3, case := 'positive correlations']

# Plots
ggplot(dt_gamma, aes(x = value, col = type, fill = type, alpha = alpha)) +
  #geom_histogram(position = 'identity', alpha = 0.1) +
  geom_density() +
  facet_wrap(~case, scales = 'free') +
  #facet_grid(row = vars(name), col = vars(id_case), scales= 'free') +
  scale_color_colorblind(labels = c("Original\nsample",
                                    "Multivariate\nsample",
                                    "Univariate\nsample")) +
  scale_fill_colorblind(labels = c("Original\nsample",
                                   "Multivariate\nsample",
                                   "Univariate\nsample")) +
  scale_alpha_continuous(limits = c(0,1), guide = 'none') +
  labs(fill = 'sample type', col = 'sample type') +
  xlab("Total steel production (sum of three disaggregate sectors)") +
  ylab("Probability density") +
  theme(legend.position = 'bottom')


# save plot
ggsave(filename = './figures_overleaf/figures/figure4.pdf',
       width = 7, height = 3)



# _b) Lognormal vs. gamma dist ===================================================
dt_lognorm <-  dt2[name %in% c('uv_lnorm', 'mv_lnorm')]

dt_lognorm <- dt_lognorm[!(name == 'uv_lnorm' & variable == 'sample')]
dt_lognorm[variable == 'sample', type := 'original']
dt_lognorm[variable == 'sample_new' & name == 'mv_lnorm', type := 'multivariate']
dt_lognorm[variable == 'sample_new' & name == 'uv_lnorm', type := 'univariate']
dt_lognorm[type == 'original', alpha := 0.6]
dt_lognorm[type != 'original', alpha := 0.3]
dt_lognorm[, type := factor(type, levels = c('original', 'multivariate', 'univariate'))]



dt_lognorm[id_case == 1, case := 'negative correlations']
dt_lognorm[id_case == 2, case := 'no correlations']
dt_lognorm[id_case == 3, case := 'positive correlations']

# Plots
ggplot(dt_lognorm, aes(x = value, col = type, fill = type, alpha = alpha)) +
  #geom_histogram(position = 'identity', alpha = 0.1) +
  geom_density() +
  facet_wrap(~case, scales = 'free') +
  #facet_grid(row = vars(name), col = vars(id_case), scales= 'free') +
  scale_color_colorblind(labels = c("Original\nsample",
                                    "Multivariate\nsample",
                                    "Univariate\nsample")) +
  scale_fill_colorblind(labels = c("Original\nsample",
                                   "Multivariate\nsample",
                                   "Univariate\nsample")) +
  scale_alpha_continuous(limits = c(0,1), guide = 'none') +
  labs(fill = 'sample type', col = 'sample type') +
  theme(legend.position = 'bottom')

ggsave(filename = './figures_overleaf/figures/figure4_lnorm.pdf',
       width = 7, height = 3)




# the end =====================================================================

