#'
#'
#'
#' @author Simon Schulte
#' Date: 2024-02-12 16:12:27.318764
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
library(ggthemes)

############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())

source("functions.R")
############################################################################## #
##### functions #################################################################
############################################################################## #


############################################################################## #
##### load data #############################################################
############################################################################## #
N <- 1E5
#alpha <- as.numeric(rdirichlet(1, rep(1,5)))
alpha <- c(0.1, 0.3, 0.6)
K <- length(alpha)


seq <- seq(0.05, 1, length.out = 12)
betas <- lapply(seq,
                function(x) alpha * x)
names(betas) <- as.character(round(seq, 2))
# betas <- list(
#   alpha * 0.01,
#   alpha * 0.1,
#   alpha * 1
# )
sample_dt <- as.data.table(expand.grid(alpha = list(alpha), beta = betas))
sample_dt[, beta_cv := pmap_dbl(list(x = alpha, y = beta), function(x,y) {
  (y / x)[1]
})]
sample_dt[, case := as.factor(1:.N)]
sample_dt[, sample := pmap(list(alpha = alpha, beta = beta), function(alpha, beta) {
  rdirg(N, alpha, beta)
})]

sample_dt[, sample := pmap(list(x = sample), function(x) {
  x <- convert_sample_to_dt(x)
  setkey(x, y)
  return(x)
})]

sample_dt2 <- sample_dt[, rbindlist(sample), by = .(case, beta_cv)]
sample_dt2[, y := as.factor(y)]
sample_dt2 <- sample_dt2[, list(value = unlist(sample, FALSE)), by = .(case, y, beta_cv)]

par_values <- sample_dt[, list(beta = unlist(beta, FALSE),
                               alpha = unlist(alpha, FALSE)),
                        by = case]
par_values[, y := as.factor(rep(1:K, length(betas)))]


ggplot(sample_dt2, aes(x = value, fill = y)) +
  #geom_histogram(position = 'identity', alpha = 0.6) +
  geom_density(alpha = 0.6, col = NA) +
  # geom_vline(data = par_values, aes(xintercept = alpha, col = y)) +
  # geom_vline(data = par_values, aes(xintercept = alpha-beta, col = y),
  #            linetype = 'dashed') +
  # geom_vline(data = par_values, aes(xintercept = alpha+beta, col = y),
  #            linetype = 'dashed') +
  facet_wrap(~as.factor(beta_cv), scales = 'free_y') +
  scale_fill_colorblind() +
  scale_color_colorblind()


sample_dt3 <- sample_dt2[, list(
  q25 = quantile(value, 0.025),
  q33 = quantile(value, 0.33),
  q67 = quantile(value, 0.67),
  q975 = quantile(value, 0.975)
), by = .(case, y)]

sample_dt3 <- merge(sample_dt3, par_values, by = c('case', 'y'))
ggplot(sample_dt3, aes(x = y)) +
  geom_errorbar(aes(ymin = q25, ymax = q975), col = 1, width = 0.3) +
  geom_errorbar(aes(ymin = q33, ymax = q67), col = 2, width = 0.3) +
  geom_errorbar(aes(ymin = alpha - beta, ymax = alpha + beta), col = 3, width = 0.3) +
  facet_wrap(~case)


sample_cv <- sample_dt2[, list(cv = sd(value) / mean(value)), by=.(case, y)]
ggplot(sample_cv, aes(x = case, y = cv, col = y)) +
  geom_point() +
  scale_color_colorblind()



############################################################################## #
##### save results #############################################################
############################################################################## #


# THE END ---------------------------------------------------------------------
