#'
#'
#'
#' @author Simon Schulte
#' Date: 2024-01-09 13:15:07.471709
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
library(Ternary)
library(cowplot)
library(corrplot)
library(patchwork)
library(grid)
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
#library(MethylCapSig)

############################################################################## #
##### settings #################################################################
############################################################################## #
options("datatable.print.class" = TRUE)
theme_set(theme_bw())
my_scale_fill <-scale_fill_colorblind()
my_cols <- (colorblind_pal()(8))

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
source('functions.R')
source('src/functions_dirichlet.R')

############################################################################## #
##### set parameters #############################################################
############################################################################## #
N <- 10000


pars_agg <- as.data.table(tribble(
  ~dist2, ~mean, ~sd, ~min, ~max,
  'unif', NA, NA, 0, 20,
  'exp', 10, NA, 0, Inf,
  'norm', 10, 8, -Inf, Inf,
  'truncnorm', 10, 8, 0, Inf
))



pars_shares <- as.data.table(tribble(
  ~dist_disagg2, ~shares, ~sds,
  'Dir(1,1,...)', c(1/3, 1/3, 1/3), NA,
  'Dir_MaxEnt', c(0.1, 0.3, 0.6), NA,
  'Dir_Generalised', c(0.1, 0.3, 0.6),  c(0.01, 0.03, 0.06),
  'Dir_Nested', c(0.1, 0.3, 0.6),  c(NA, NA, 0.01)
))



############################################################################## #
##### plot #############################################################
############################################################################## #

# standard plot
plot1_full(pars_agg, pars_shares, colour_headers = 'grey90')

ggsave(filename = './figures_overleaf/figures/figure3.pdf',
       width = 8, height = 8)



############################################################################## #
##### save results #############################################################
############################################################################## #


# THE END ---------------------------------------------------------------------
