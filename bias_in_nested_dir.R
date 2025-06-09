#'
#'
#'
#' @author Simon Schulte
#' Date: 2025-05-25 16:38:58.51824
#'
#' Content:
#'


############################################################################## #
##### load packages ############################################################
############################################################################## #

library(data.table)
library(tidyverse)
library(MaxentDisaggregation)
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


############################################################################## #
##### functions #################################################################
############################################################################## #


############################################################################## #
##### load data #############################################################
############################################################################## #
sample_shares <- rshares(1E5,
                         c(a=0.2,b=0.3,c=NA, d= NA),
                         c(a=1E-5,b=0.05, c=NA, d = NA))


plot_samples(sample_shares) +
  facet_wrap(~variable)
colMeans(sample_shares)
round(colSds(sample_shares), 3)

set.seed(1)
shares <- as.numeric(rdirichlet_uniform(1, 100))
sds <- shares * as.numeric(runif(100, min = 0, max = 1.3))

sds[1:20] <- NA
sample_shares <- rshares(1E5, shares, sds)
# plot_samples(sample_shares) +
#   facet_wrap(~variable) +
#   theme(legend.position = 'none')

data.table(sample_sd = colSds(sample_shares),
           provided_sd = sds) %>%
  na.omit %>%
  ggplot(aes(x = provided_sd, y = sample_sd)) +
  geom_point()+
  geom_abline(slope = 1, intercept = 0, col = 'red')

data.table(sample_mean = colMeans(sample_shares),
           provided_mean = shares,
           provided_cv = sds / shares) %>%
  na.omit %>%
  ggplot(aes(x = provided_mean, y = sample_mean)) +
  geom_point(alpha = 0.3, aes(size = provided_cv))+
  geom_abline(slope = 1, intercept = 0, col = 'red')

data.table(sample_mean = colMeans(sample_shares),
           provided_mean = shares,
           provided_cv = sds / shares) %>%
  na.omit %>%
  ggplot(aes(x = provided_cv, y = (sample_mean - provided_mean) / provided_mean)) +
  geom_point(alpha = 0.3, aes(size = provided_cv))+
  geom_hline(yintercept = 0, col = 'red')

data.table(sample_sd = colSds(sample_shares),
           provided_sd = sds,
           provided_mean = shares) %>%
  na.omit %>%
  ggplot(aes(x = provided_sd, y = (sample_sd - provided_sd) / provided_sd)) +
  geom_point(alpha = 0.3)+
  geom_hline(yintercept = 0, col = 'red')




############################################################################## #
##### save results #############################################################
############################################################################## #


# THE END ---------------------------------------------------------------------
