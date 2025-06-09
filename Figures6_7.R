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
# library(MethylCapSig)
library(matrixStats)
library(MaxentDisaggregation)
library(GGally)
library(patchwork)
library(rlang)

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


#' Draw N samples from k independent log-normal distributions
#'
#' @param N        integer, number of samples (rows of the output matrix)
#' @param meanlog  numeric vector of log-means (µ₁ … µ_k)
#' @param sdlog    numeric vector of log-SDs   (σ₁ … σ_k)
#' @return         N × k matrix; column j contains rlnorm(N, meanlog[j], sdlog[j])
sample_lognorm <- function(N, meanlog, sdlog)
{
  stopifnot(
    length(meanlog) == length(sdlog),
    N > 0, is.finite(N)
  )

  # draw each column independently
  mat <- sapply(seq_along(meanlog), function(j)
    rlnorm(N, meanlog = meanlog[j], sdlog = sdlog[j])
  )

  # make sure it is a matrix even if k == 1
  mat <- as.matrix(mat)
  colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  mat
}
my_lower <- function(data, mapping, ...) {
  # filled density
  ggally_density(
    data   = data,
    alpha = 0.6,
    mapping = modifyList(mapping,     # keep x & y, but add the fill aesthetic
                         aes(fill = after_stat(level))),
    ...
  ) +
    ## add contour lines on top
  # stat_density_2d(                     # only the lines
  #   data     = data,
  #   mapping  = mapping,
  #   #contour = FALSE,
  #   # same x and y
  #   colour   = "black",
  #   linewidth = .3,                    # thin lines look nicer in small panels
  #   ...) +
   ggally_points(data = data, mapping = mapping, ...) +
    scale_fill_viridis_c()
}

## helper for lower-triangle panels ----------------------------------------
lower_pct <- function(data, mapping, ...) {
  # column names that ggpairs mapped to x and y
  x_col <- as_name(mapping$x)
  y_col <- as_name(mapping$y)

  # percentage of points above the 45° line
  pct <- mean(data[[y_col]] > data[[x_col]], na.rm = TRUE) * 100

  ggplot(data = data, mapping = mapping) +
    geom_point(size = 0.1, alpha = 0.2, colour = "grey30") +        # raw points
    geom_abline(slope = 1, intercept = 0, colour = "red") +         # y = x line
    annotate(
      "text",
      x = Inf, y =  Inf,                                           # top-left corner
      hjust = 1.2, vjust =  1.2,
      label = sprintf("P(Y>X) = %.0f%%", pct),                               # formatted %
      size  = 3.2
    ) +
    theme_bw()
}

## build the matrix --------------------------------------------------------
# ggpairs(
#   as.data.table(param_dt$sample[[1]][1:1000, ]) |>
#     setnames(c("HBEV", "BEV", "ICE")),
#
#   upper = list(continuous = function(...)
#     my_fn(..., col = "black",
#           title = "Corr", digits = 2, stars = FALSE)),
#
#   lower = list(continuous = lower_pct),   # <-- use the new function here
#
#   diag  = list(
#     continuous = function(data, mapping, ...)
#       diag_hist_colour(data, mapping,
#                        palette = my_cols, alpha = 0.30, ...)),
#
#   title = "Disaggregates"
# )

# ggpairs(as.data.table(param_dt$sample[[1]][1:1000,]) %>%
#           setnames(c('HBEV', 'BEV', 'ICE')),
#         upper = list(continuous = function(...) my_fn(..., col = 'black',
#                                                       title = 'Corr', digits = 2,
#                                                       stars = FALSE)),
#         #lower = list(continuous = wrap('points', )),
#         #lower = list(continuous = my_lower),
#         lower = list(continuous = function(...) ggally_points(..., size = 0.3,
#                                                               alpha = 0.4,
#                                                               col = 'grey30',
#                                                               shape = 16)+
#                        geom_abline(slope = 1, intercept = 0, col = 'red') +
#                        annotate(geom = 'text', x = 25, y = 25,
#                                 label = "plot here the percentage of points above the red line") +
#                        theme_bw()),
#         diag = list(
#           continuous = function(data, mapping, ...)
#             diag_hist_colour(data, mapping,
#                              palette = my_cols,   # remove this line to use default hues
#                              alpha = 0.3, ...)
#         ),
#         title = 'Disaggregates'
# )
# as.data.table(param_dt$sample[[1]][1:1000,]) %>%
#   setnames(c('HBEV', 'BEV', 'ICE')) %>%
#   .[BEV > ICE, .N / N]

############################################################################## #
##### load data #############################################################
############################################################################## #

# 1. Set parameters ===============================================================
N <- 1E5
bins <- 50

param_dt <- as.data.table(tribble(
  ~id_case, ~mean_0, ~sd_0, ~shares, ~sds,
  1, 100, 5, c(0.1, 0.3, 0.6), NULL,
  2, 100, 58, c(NA, NA, 0.6), c(NA, NA, 0.03),
  #2, 100, 10, c(0.1, 0.3, 0.6), c(0.1, 0.3, 0.6)*0.3,
  3, 100, 30, c(0.1, 0.3, 0.6), c(0.1, 0.3, 0.6)*0.005
))

param_dt <- as.data.table(tribble(
  ~id_case, ~mean_0, ~sd_0, ~shares, ~sds,
  1, 100, 5, c(0.2, 0.35, 0.45), NULL,
  2, 100, 58, c(NA, NA, 0.6), c(NA, NA, 0.03),
  #2, 100, 10, c(0.1, 0.3, 0.6), c(0.1, 0.3, 0.6)*0.3,
  3, 100, 30, c(0.2, 0.38, 0.42), c(0.2, 0.38, 0.42)*0.005
))
# 2. Sample ===============================================================

# _a) sample agg, shares and disaggs ==========================================
param_dt[, sample := pmap(list(mean_0 = mean_0, sd_0 = sd_0,
                               shares = shares, sds = sds, min = 0, n= N,
                               log = TRUE),
                          rdisagg)]

param_dt[, sample_agg := pmap(list(mean = mean_0, sd = sd_0, min = 0, n= N,
                                   log = TRUE),
                              ragg)]

param_dt[, sample_shares := pmap(list(shares = shares, sds = sds, n= N,
                                      log = TRUE),
                                 rshares)]

# _b) correlaionts ============================================================

param_dt[, cor := pmap(list(x = sample), cordt)]

param_dt[, .(id_case, mean_0, sd_0, shares, sds)]
param_dt$cor


# _c) reshape dt ================================================================
param_dt2 <- melt(param_dt, id.vars = c('id_case'),
                  measure.vars = c('sample', 'sample_agg', 'sample_shares', 'cor'))

samples_agg <- param_dt$sample_agg
samples_agg <- lapply(samples_agg, function(x) data.table(value = x))
names(samples_agg) <- param_dt$id_case
samples_agg <- rbindlist(samples_agg, idcol = 'id_case')

samples_shares <- param_dt$sample_shares
samples_shares <- lapply(samples_shares, function(x) {
  x <- reshape2::melt(x)
  x <- as.data.table(x)
  x
})
names(samples_shares) <- param_dt$id_case
samples_shares <- rbindlist(samples_shares, idcol = 'id_case')

# _d) Resample disaggregates from univariate and multivariate dists ============
# __ univariate lognorm =======================================================
param_dt[, meanlog := lapply(sample, function(x) x %>% log %>% colMeans)]
param_dt[, sdlog := lapply(sample, function(x) x %>% log %>% colSds)]
param_dt[, sample_nocorr := Map(sample_lognorm, meanlog = meanlog, sdlog = sdlog,
                                N = N)]
sample_nocorr <- param_dt$sample_nocorr
sample_nocorr <- lapply(sample_nocorr, function(x) {
  x <- reshape2::melt(x)
  x <- as.data.table(x)
  x
})
names(sample_nocorr) <- param_dt$id_case
sample_nocorr <- rbindlist(sample_nocorr, idcol = 'id_case')

sample_nocorr_agg <- sample_nocorr[, .(value = sum(value)), by = .(id_case, Var1)]

# __ univariate gamma =========================================================
param_dt[, mean := lapply(sample, function(x) x %>%  colMeans)]
param_dt[, sd := lapply(sample, function(x) x %>%  colSds)]
param_dt[, sample_gamma := Map(rgamma2_vec, mean = mean, sd = sd,
                               n = N)]
sample_gamma <- param_dt$sample_gamma
sample_gamma <- lapply(sample_gamma, function(x) {
  x <- reshape2::melt(x)
  x <- as.data.table(x)
  x
})
names(sample_gamma) <- param_dt$id_case
sample_gamma <- rbindlist(sample_gamma, idcol = 'id_case')
sample_gamma_agg <- sample_gamma[, .(value = sum(value)), by = .(id_case, Var1)]

# __multivariate Gamma/lognorm =========================================================

dists_dt <- as.data.table(tribble(
  ~id_dist, ~ dist, ~ type, ~name,
  #  1, 'lnorm', 'univariate2', 'uv_lnorm',
  2, 'lnorm', 'rlnorm.rplus', 'mv_lnorm',
  #3, 'gamma', NULL, 'uv_gamma',
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

# Unnest =================================
dt2 <- dt[, .(id_dist, id_case, name, sample_agg, sample_agg_new)] %>%
  .[, list(sample = unlist(sample_agg, recursive = FALSE),
           sample_new = unlist(sample_agg_new, recursive = FALSE)),
    by = .(id_dist, id_case, name)]

dt2 <- melt(dt2, id.vars = c('id_dist', 'id_case', 'name'))

sample_mvgamma_agg <- dt2[name == 'mv_gamma' & variable == 'sample_new',
                          .(id_case, value)]
sample_mvlognorm_agg <- dt2[name == 'mv_lnorm' & variable == 'sample_new',
                            .(id_case, value)]




# _ Combine samples
samples_agg_combined <- rbindlist(
  list(
    original = samples_agg,
    lognorm = sample_nocorr_agg,
    gamma = sample_gamma_agg,
    mvgamma = sample_mvgamma_agg,
    mvlognorm = sample_mvlognorm_agg
  ),
  idcol = 'type', use.names = TRUE, fill = TRUE
)

samples_agg_combined[, type := factor(type, levels = c('original', 'lognorm', 'gamma',
                                                       'mvlognorm', 'mvgamma'))]

# 3. Plots  ============================================================

# Figure 5 =------------------------------------------------------------
# case 1
p_agg <- lapply(1:nrow(param_dt), function(i) {
  ggplot(samples_agg[id_case == i], aes(x = value)) +
    geom_histogram(alpha = 0.3, col = 'black', fill = 'black', bins = bins) +
    ggtitle(paste0('Aggregate: ',
                   fmt_scalar_for_title(param_dt[id_case == i, mean_0][[1]], 'm₀'),
                   ', ',
                   fmt_scalar_for_title(param_dt[id_case == i, sd_0][[1]], 's₀')))

})
p_shares <- lapply(1:nrow(param_dt), function(i) {
  ggplot(samples_shares[id_case == i],
         aes(x = value, col = as.factor(Var2), fill = as.factor(Var2))) +
    geom_histogram(alpha = 0.3, position = 'identity', bins = bins) +
    scale_color_colorblind7() +
    scale_fill_colorblind7() +
    theme(legend.position = 'none') +
    ggtitle(paste0('Shares: ',
                   fmt_vec_for_title(param_dt[id_case == i, shares][[1]], '𝐦'),
                   ',\n\t       ',
                   fmt_vec_for_title(param_dt[id_case == i, sds][[1]], '𝐬', digits = 4)))



})



my_cols <- colorblind_pal()(8)[2:4]
# names(my_cols) <-  as.character(1:3)
names(my_cols) <-  c('HBEV', 'BEV', 'ICE')
plist <- lapply(param_dt$sample, function(x) {
  #sample_disagg[[1]] %>% as.data.table
  GGally::ggmatrix_gtable(
    ggpairs(as.data.table(x) %>% setnames(c('HBEV', 'BEV', 'ICE')),
            upper = list(continuous = function(...) my_fn(..., col = 'black',
                                                          title = 'Corr', digits = 2,
                                                          stars = FALSE)),
            #lower = list(continuous = wrap('points', )),
            #lower = list(continuous = my_lower),
            lower = list(continuous = lower_pct),   # <-- use the new function here
            #
            # lower = list(continuous = function(...) ggally_points(...,
            #                                                       size = 0.1,
            #                                                       alpha = 0.1,
            #                                                       col = 'grey30',
            #                                                       shape = 16)+
            #                theme_bw()),
            diag = list(
              continuous = function(data, mapping, ...)
                diag_hist_colour(data, mapping,
                                 palette = my_cols,   # remove this line to use default hues
                                 alpha = 0.3, ...)
            ),
            title = 'Disaggregates'
    )
  )
})
plot(plist[[1]])



p_full <- list()
for (i in 1:nrow(param_dt)) {
  p_full[[i]] <- wrap_elements(
    ((p_agg[[i]] / p_shares[[i]]) | plist[[i]] ) +
      plot_layout(widths = c(1.5, 2)) +
      plot_annotation(title = paste0('Scenario ', i),
                      theme = theme(plot.title = element_text(size = 16)))
  )

}

p_full[[1]]

wrap_plots(p_full, ncol = 1)

ggsave(filename = './figures_overleaf/figures/figure_example1.png',
       width = 7, height = 9, scale  = 1.3)

ggsave(filename = './figures_overleaf/figures/figure_example1.pdf',
       width = 7, height = 9, scale  = 1.3, device = cairo_pdf)




# Figure 6 ====================================================================
samples_agg_combined[id_case == '1', Scenario := "1) Negative Correlations"]
samples_agg_combined[id_case == '2', Scenario := "2) Mixed Correlations"]
samples_agg_combined[id_case == '3', Scenario := "3) Positive Correlations"]

emission_factor_steel <- 2.5
samples_agg_combined[, value := value * emission_factor_steel]

summary_stats <- samples_agg_combined[type %in% c('original', 'gamma', 'mvgamma'),
                     .(Scenario,
                       value = value,
                       type = type)] %>%
  .[, .(
    Mean = mean(value),
    "2.5th perc" = quantile(value, 0.025),
    "97.5th perc" = quantile(value, 0.975)
  ), by = .(type, Scenario)]

summary_stats <- melt(summary_stats, id.vars = c('type', 'Scenario'),
                      variable.name = 'Information shared\nby Model 1')

# p_nocorr <-   ggplot(samples_agg_combined[type %in% c('original', 'gamma', 'mvgamma'),
#                                           .(Scenario = id_case,
#                                             value = value,
#                                             type = type)],
#                      aes(x = value, col = type, fill = type)) +
#   geom_histogram(alpha = 0.5, bins = bins, linewidth = 0.15,
#                  position = 'identity') +
#   geom_vline(data = summary_stats, aes(xintercept = value,
#                                        linetype = `Information shared\nby Model 1`,
#                                        col = type)) +
#   scale_color_colorblind() +
#   scale_fill_colorblind() +
#   labs(col = 'Sample\ntype', fill = 'Sample\ntype') +
#   facet_wrap(~Scenario, labeller = 'label_both', scales = 'free_x') +
#   theme(legend.position = 'bottom')
#
# p_nocorr

p_nocorr1 <-   ggplot(samples_agg_combined[type %in% c('original', 'gamma', 'mvgamma'),
                                          .(Scenario,
                                            value = value,
                                            type = type)],
                     aes(x = value, col = type, fill = type)) +
  geom_histogram(alpha = 0.5, bins = bins, linewidth = 0.15,
                 position = 'identity') +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  labs(col = 'Information shared\nby Model 1', fill = 'Information shared\nby Model 1') +
  coord_cartesian(xlim = c(0,750)) +
  facet_wrap(~Scenario, labeller = 'label_both', scales = 'free_x') +
  theme(legend.position = 'none', axis.title.x = element_blank())


p_nocorr2 <-   ggplot(samples_agg_combined[type %in% c('original', 'gamma', 'mvgamma'),
                                           .(Scenario,
                                             value = value,
                                             type = type)],
                      aes(x = value,y = type, col = type, fill = type)) +
  # geom_boxplot(alpha = 0.3, outlier.alpha = 0.3, outlier.shape = 16,
  #              outlier.size = 0.8) +
  geom_point(alpha = 0, size = 0) +
  stat_summary(fun.data =  median.quartile, geom = 'crossbar', width = .4,
               fill = 'white') +
  # geom_jitter(shape=16, alpha = 0.2, col = 'grey50', size = 0.6,
  #             position=position_jitter(0.1, seed = 1)) +
  stat_summary(fun.data =  percentile, geom = 'errorbar', linewidth = .4,
               size = 0.2, width = 0.2) +
  stat_summary(fun.data =  median.quartile,
               geom = 'crossbar', width = .4, size = 0.7,
               alpha=0.4) +
  scale_color_colorblind(
    name   = "Information shared\nby Model 1",
    breaks = c("original", "gamma", "mvgamma"),
    labels = c("Full sample", "Mean & Standard Deviation",
               "Mean + Covariance Matrix")
  ) +
  scale_fill_colorblind(
    name   = "Information shared\nby Model 1",
    breaks = c("original", "gamma", "mvgamma"),
    labels = c("Full sample", "Mean & Standard Deviation",
               "Mean + Covariance Matrix")
  ) +
  labs(col = 'Information shared\nby Model 1', fill = 'Information shared\nby Model 1') +
  facet_wrap(~Scenario, labeller = 'label_both', scales = 'free_x') +
  coord_cartesian(xlim = c(0,750)) +
  #theme_minimal() +
  #theme_void() +
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank() ,
        strip.background = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ))
# p_nocorr2

p_nocorr1 / p_nocorr2 +
  plot_layout(heights = c(3,1)) +
  plot_annotation(tag_levels = 'a')

ggsave(filename = './figures_overleaf/figures/figure_example2.pdf',
       width = 7, height = 5, scale  = 1)


ggsave(filename = './figures_overleaf/figures/figure_example2.png',
       width = 7, height = 5, scale  = 1)




# the end ====================================================================
samples_agg_combined[, .(
  mean = mean(value) %>% round(2),
  sd = sd(value) %>% round(2),
  cv = (sd(value) / mean(value)) %>% round(2),
  q025 = quantile(value, 0.025) %>% round(2),
  q5 = quantile(value, 0.5) %>% round(2),
  q975 = quantile(value, 0.975) %>% round(2)

), by = .(type, id_case)]

# the end =====================================================================

