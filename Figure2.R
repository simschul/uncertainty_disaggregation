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
#library(ggforce)
library(gtools)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
#library(faux)
# library(compositions)
library(disaggR)
library(ggthemes)
library(TruncatedNormal)
# library(moments)
# library(mvLognCorrEst)
#library(MethylCapSig)
library(matrixStats)
library(MaxentDisaggregation)
library(GGally)
library(ggthemes)
library(matrixStats)
library(grid)
library(patchwork)
library(rlang)      # for the %||% operator


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
source('src/functions_plot.R')

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


make_vert_label <- function(label, base_theme = ggplot2::theme_get(),
                            angle = 90, hjust = 0.5, vjust = 0.5) {

  ## 1. resolve the *effective* plot.title element

  ##    calc_element() turns any rel() units into points and
  ##    merges theme defaults with your overrides
  title_el <- ggplot2:::calc_element("plot.title", base_theme)

  ## 2. translate that element into a grid::gpar
  gp <- gpar(
    fontsize   = title_el$size,
    fontface   = title_el$face   %||% "plain",
    fontfamily = title_el$family %||% "",
    col        = title_el$colour %||% "black",
    lineheight = title_el$lineheight %||% 1
  )

  ## 3. build a rotated text grob and wrap it so patchwork sees one leaf
  wrap_elements(
    full = textGrob(label, gp = gp, rot = angle,
                    hjust = hjust, vjust = vjust)
  )
}


myplot <- function(sample, log = FALSE) {
  sample_dt <- as.data.table(sample)
  setnames(sample_dt, c('A', 'B'))
  sample_dt[A > B, .N/N]

  my_cols <- colorblind_pal()(8)[2:3]
  names(my_cols) <- c("A", 'B')


  sample_long <- melt(sample_dt)
  sample_summary <- sample_long[, .(
    Mean = mean(value),
    "2.5th perc" = quantile(value, 0.025),
    '97.5th perc' = quantile(value, 0.975)
  ), by = variable]
  sample_summary <- melt(sample_summary ,id.vars = c("variable"),
                         variable.name = "Sample\nstatistics")
  probs <- data.table(
    prop = round(sample_dt[B > A, .N/N], 2),
    x = median(sample_long$value),
    x2 = max(sample_long$value) / 2
  )

  if (log == TRUE) {
    probs[, x := log10(x)]
  }

  plist <- list()
  plist[[2]] <- ggplot(sample_dt, aes(x = A, y = B)) +
    geom_abline(intercept = 0, slope = 1, col="red") +
    geom_point(shape = 16, size = 0.5, alpha = 0.6) +
    annotate(geom = 'text',
             x = probs$x, y = Inf,
             label = paste0('P(B>A) = ', probs$prop),
             vjust = 2.6)

  if (isTRUE(log)) {
    plist[[2]] <- plist[[2]] + scale_x_log10() + scale_y_log10()
  }

  plist[[1]]<- ggplot(sample_long, aes(x = value, col = variable, fill = variable)) +
    geom_histogram(alpha = 0.3, position = 'identity', bins = 60) +
    # geom_vline(data = sample_summary,
    #            aes(xintercept = mean, col = variable),
    #            linetype="dotted") +
    geom_vline(data = sample_summary,
               aes(xintercept = value, col = variable,
                   linetype = `Sample\nstatistics`)) +


    # annotate(geom = 'text',
    #          x = probs$x2, y = Inf,
    #          label = paste0('P(B>A) = ', probs$prop),
    #          vjust = 2.6) +
    ylab("Count") +
    labs(fill = 'Variable', col = 'Variable') +
    scale_linetype_manual(values=c("dotted", "dashed", "twodash")) +
    scale_fill_colorblind7() +
    scale_color_colorblind7()

  return(plist)

}


myplot2 <- function(dt_gamma) {
  dt_gamma_summary <- dt_gamma[, .(
    Mean = mean(value),
    "2.5th perc" = quantile(value, 0.025),
    '97.5th perc' = quantile(value, 0.975)
  ), by = .(case, name, type)]
  dt_gamma_summary <- melt(dt_gamma_summary,
                           id.vars = c('case', 'name', 'type'))
  # dt_gamma_summary[, variable := factor(variable,
  #                                       levels = c('', 'Mean',
  #                                                  "2.5th perc",
  #                                                  '97.5th perc'))]

  ggplot(dt_gamma[type != 'multivariate' & case != 'no correlations'],
         aes(x = value, col = type, fill = type, alpha = alpha)) +
    geom_histogram(aes(col = type, fill = type),
                   position = 'identity', alpha = 0.3, bins = 60) +
    geom_vline(data = dt_gamma_summary[type != 'multivariate'
                                       & case != 'no correlations'] ,
               aes(xintercept = value, col = type, linetype = variable)) +
    scale_color_colorblind2(labels = c("Original\nsample",
                                       "Ignoring\ncorrelations")) +
    scale_fill_colorblind2(labels = c("Original\nsample",
                                      "Ignoring\ncorrelations")) +
    scale_alpha_continuous(limits = c(0,1), guide = 'none') +
    scale_linetype_manual(values=c("dotted", "dashed", "twodash")) +
    labs(fill = 'Sample type', col = 'Sample type',
         linetype = "Sample\nstatistics") +
    xlab("A + B") +
    ylab("Count") +
    theme(legend.position = 'bottom')
}

############################################################################## #
##### Problem 1: decision making #############################################################
############################################################################## #

N <- 1E5
# product comparison with correlations
sample_agg <- ragg(N, mean = 100, sd = 80, log = FALSE, min = 0)
sample_shares <- rshares(N, shares = c(0.45, 0.55), sds = c(1E-3, 1E-3))
sample <- sample_shares * sample_agg
colSds(sample) / colMeans(sample)

p_pos <- myplot(sample, log = TRUE)

# ==============================================================================
sample_agg <- ragg(N, mean = 100, sd = 5, log = TRUE)
sample_shares <- rshares(N, shares = c(0.45, 0.55))
sample <- sample_shares * sample_agg
colSds(sample) / colMeans(sample)


p_neg <- myplot(sample)



############################################################################## #
##### Problem 2: model 2 uncertainties #############################################################
############################################################################## #


# 1. Set parameters ===============================================================
N <- 1E5
param_dt <- as.data.table(tribble(
  ~id_case, ~mean_0, ~sd_0, ~shares, ~sds,
  1, 100, 5, c(0.4, 0.6), NULL,
  2, 100, 58, c(0.5, 0.5), NULL,
  3, 100, 30, c(0.49, 0.51), c(1E-3, 1E-3)
))

# 2. Sample ===============================================================

param_dt[, sample := pmap(list(mean_0 = mean_0, sd_0 = sd_0,
                               shares = shares, sds = sds, min = 0, n= N,
                               log = TRUE),
                          rdisagg)]

param_dt$sample[[1]] %>% rowSums()
param_dt[, cor := pmap(list(x = sample), cordt)]

# 3. Resample ===============================================================
dists_dt <- as.data.table(tribble(
  ~id_dist, ~ dist, ~ type, ~name,
  #  1, 'lnorm', 'univariate2', 'uv_lnorm',
  #  2, 'lnorm', 'mvlognormal', 'mv_lnorm',
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

plist2 <- list()
# Plots
plist2[['neg']] <- myplot2(dt_gamma[case == 'negative correlations'])
plist2[['pos']] <- myplot2(dt_gamma[case == 'positive correlations'])



# assemble plots ===============================================================
p_pos[[1]] <- p_pos[[1]] +
  ggtitle("Marginal histograms of A and B") +
  theme(axis.title.x = element_blank())
p_pos[[2]] <- p_pos[[2]] +
  ggtitle("Correlation between A and B") +
  theme(axis.title.x = element_blank())
plist2$pos <- plist2$pos +
  ggtitle("Implication of ignoring cor-\nrelations on the sum of A and B") +
  theme(axis.title.x = element_blank())


# 2. create vertical label grobs

lab_pos <- make_vert_label("Positive correlation")
lab_neg <- make_vert_label("    Negative correlation")

plots_left  <- (p_pos[[1]] + p_pos[[2]] + plist2$pos)  /
  (p_neg[[1]] + p_neg[[2]] + plist2$neg)

labels_right <- lab_pos / lab_neg



(labels_right | plots_left) +
  # 2 columns at top level
  plot_layout(widths = c(0.05, 1), guides = "collect") &
  theme(legend.position = "bottom")


ggsave(filename = './figures_overleaf/figures/figure1b.pdf',
       width = 7, height = 4, scale = 1.6)
ggsave(filename = './figures_overleaf/figures/figure1b.png',
       width = 7, height = 4, scale = 1.6)


# save plot
