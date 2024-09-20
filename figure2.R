#'
#'
#'
#' @author Simon Schulte
#' Date: 2024-04-11 16:10:37.185081
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


############################################################################## #
##### functions #################################################################
############################################################################## #
source('functions.R')






############################################################################## #
##### load data #############################################################
############################################################################## #
N <- 1E5
shares <- c(0.1, 0.3, 0.6)
gammas <- c(1, 6.36, 30)

sample <- lapply(gammas, function(x) rdir(N, shares, gamma = x))
dt <- lapply(sample, melt_matrix) %>%
  setNames(gammas) %>%
  rbindlist(idcol = 'gamma')

dt[, gamma := as.factor(gamma)]
dt[, gamma := paste0('gamma == ', gamma)]
dt[, gamma := factor(gamma, levels = c('gamma == 1', 'gamma == 6.36', 'gamma == 30'))]
#dt[gamma == 'gamma == 6.36', gamma := 'gamma == 6.36 (MaxEnt)']
dt[, y := as.character(y)]
ggplot(dt, aes(x = value, col = y, fill = y)) +
  geom_histogram(alpha = 0.3, position = 'identity') +
  facet_wrap(~gamma, labeller = label_parsed) +
  labs(fill = 'variable', col = 'variable') +
  scale_fill_colorblind7(labels = c('1' = bquote(~Y[1]),
                                    '2' = bquote(~Y[2]),
                                    '3' = bquote(~Y[3]))) +
  scale_color_colorblind7(labels = c('1' = bquote(~Y[1]),
                                     '2' = bquote(~Y[2]),
                                     '3' = bquote(~Y[3])))

ggsave(filename = './figures_overleaf/figures/figure2.pdf',
       width = 6, height = 3)
############################################################################## #
##### save results #############################################################
############################################################################## #






# THE END ---------------------------------------------------------------------
