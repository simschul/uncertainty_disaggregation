#'
#'
#'
#' @author Simon Schulte
#' Date: 2024-02-08 11:34:59.678021
#'
#' Content:
#'


############################################################################## #
##### load packages ############################################################
############################################################################## #

library(data.table)
library(tidyverse)
#library(units)
library(ggforce)
library(qs)
#library(mRio)
#library(Matrix)
#library(furrr)
library(ggthemes)
library(disaggR)
library(ggridges)
library(ggrepel)
library(scales)

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

path2plots <- './figures_overleaf/figures'
path2output <- "/home/simon/Documents/Projects/Uncertainty_Extensions/code/uncertainty_GHG_accounts/intermediate_results/V1.3/"
path2exiobaseIOT <- "/home/simon/Documents/PhD_PROSET/data/EXIOBASE3/V3.8.2/IOT_2015_ixi"


############################################################################## #
##### functions #################################################################
############################################################################## #
source('functions.R')
source('./src/functions_plot.R')

geom_violin_box_jitter <- list(
  stat_summary(fun.data =  median.quartile, geom = 'crossbar', width = .2,
               fill = 'white'),
  geom_jitter(aes(size = mean %>% log10, alpha = mean %>% log10), shape=16,
              position=position_jitter(0.1, seed = 1)),
  stat_summary(fun.data =  percentile, geom = 'errorbar', width = .2,
               size = 1, col = my_cols[3]),
  geom_violin(trim = TRUE, alpha = 0.2, fill = my_cols[7], col = NA),
  stat_summary(fun.data =  median.quartile,
               geom = 'crossbar', width = .2, size = 0.7,
               alpha=0.4, fill = my_cols[3], col = my_cols[3]),

  ylab('CV = sd / mean') ,
  xlab(''),
  labs(size = 'Share of total', alpha = 'Share of total'),
  theme(legend.position = 'bottom')
)

############################################################################## #
##### load data #############################################################
############################################################################## #
fpm_sample <- qread('./intermediate_results/german_aea_fp_multiplier.qs')
fpn_sample <- qread('./intermediate_results/german_aea_fp_national.qs')
fpm_cor_sample <- qread('./intermediate_results/german_aea_fp_multiplier_cor.qs')

fpm_resampled <- qread('./intermediate_results/german_aea_fp_multiplier_resampled.qs')
fpn_resampled <- qread('./intermediate_results/german_aea_fp_national_resampled.qs')



aea_raw <- qread('./intermediate_results/german_aea_sample2.qs')
aea_cor <- qread('./intermediate_results/german_aea_cor.qs')


indices <- qread('./intermediate_results/german_aea_indices_F.qs')
Ymat <- readRDS(file.path(path2output, 'prepare_EXIOBASE_Y.RData'))
rownames(Ymat) <- as.character(1:nrow(Ymat))
Ymat2 <- matsbyname::Matrix(Ymat)
Ymat2 <- Ymat2[indices[region == 'DE']$col, 'DE', drop = FALSE]

############################################################################## #
##### Prepare data #############################################################
############################################################################## #
# omit case 5
fpm_sample <- fpm_sample[!(case %in% c(1,5))][, case := case - 1]
fpn_sample <- fpn_sample[!(case %in% c(1,5))][, case := case - 1]
fpm_cor_sample <- fpm_cor_sample[!(case %in% c(1,5))][, case := case - 1]

fpm_resampled <- fpm_resampled[!(case %in% c(1,5))][, case := case - 1]
fpn_resampled <- fpn_resampled[!(case %in% c(1,5))][, case := case - 1]

aea_raw <- aea_raw[!(case %in% c(1,5))][, case := case - 1]
aea_cor <- aea_cor[!(case %in% c(1,5))][, case := case - 1]

# 1. AEA =======================================================================

aea <- aea_raw[, rbindlist(sample), by = case]
N <- aea$sample[[1]] %>% length
aea_summary <- calculate_summary_statistics(aea)
aea_sample <- aea[, list(value = unlist(sample, FALSE)), by = .(case, y)]

aea_sample[, run_id := rep(1:N, nrow(aea))]


aea_summary <- merge(aea_summary, indices[region == 'DE',
                                          .(col = (col), CodeNr)],
                     by.x = 'y', by.y = 'CodeNr', sort = FALSE, all.x = TRUE)

#
#
#
# aea_cor <- aea_sample[, list(data = list(data.table(run_id, y, value))),
#                                 by = case]
# aea_cor[, data := pmap(list(x = data), function(x) dcast(x, run_id ~y,
#                                                                 value.var = 'value'))]
# aea_cor[, cor := pmap(list(x = data), function(x) cordt(x[, -'run_id'],
#                                                                nThreads = n_cores))]
#
# aea_cor <- aea_cor[, rbindlist(cor), by = case]

# 2. National footprint ========================================================

fpn_summary <- fpn_sample[, list(
  mean = mean(value),
  sd = sd(value),
  cv = sd(value) / mean(value),
  CI25 = quantile(value, 0.025),
  CI975 = quantile(value, 0.975),
  median = median(value)
), by = case]

# 3. Footprint multiplier =====================================================

fpm_summary <- fpm_sample[, list(
  mean = mean(value),
  sd = sd(value),
  cv = sd(value) / mean(value),
  CI25 = quantile(value, 0.025),
  CI975 = quantile(value, 0.975),
  median = median(value)
), by = .(case, j)]
fpm_summary[, mean_rel := mean / sum(mean), by = case]
fpm_summary <- merge(fpm_summary, indices[region == 'DE',
                                          .(col = as.character(col), CodeNr, sector)],
                     by.x = 'j', by.y = 'col', sort = FALSE, all.x = TRUE)

fpm_summary <- merge(fpm_summary, convert_matsbyname_to_dt(Ymat2)[, .(i, fd = value)],
                     by.x = c('j'), by.y = 'i', all.x = TRUE, sort = FALSE)

fpm_summary[, abs_mean := mean * fd]
fpm_summary[, abs_mean_rel := abs_mean / sum(abs_mean), by = case]

fpm_cv_summary <- fpm_summary[, list(median = median(cv),
                                     CI25 = quantile(cv, 0.025),
                                     CI975 = quantile(cv, 0.975)),
                              by= case]


# 4. Correlations =============================================================
#fpm_cor_sample <- fpm_cor[, rbindlist(cor), by = case]
fpm_cor_summary <- fpm_cor_sample[, list(median = median(value),
                                         CI25 = quantile(value, 0.025),
                                         CI975 = quantile(value, 0.975)),
                                  by= case]

fpm_cor_sample <- merge(fpm_cor_sample, fpm_summary[, .(j, case, cv, abs_mean,
                                                        abs_mean_rel)],
                        by.x = c('case', 'i'), by.y = c('case', 'j')) %>%
  merge(fpm_summary[, .(j, case, cv, abs_mean, abs_mean_rel)],
        by.x = c('case', 'j'), by.y = c('case', 'j'),
        suffix = c('_i', '_j'))


aea_cor2 <- merge(aea_cor, aea_summary[, .(case, y, cv, mean)],
                  by.x = c('case', 'i'), by.y = c("case", 'y'), all.x = TRUE) %>%
  merge(aea_summary[, .(case, y,cv,  mean)],
        by.x = c('case', 'j'), by.y = c("case", 'y'),
        all.x = TRUE, suffix = c('_i', '_j')) %>%
  .[, mean := pmap_dbl(list(x = mean_i, y = mean_j),
                       function(x,y) mean(c(x,y)))] %>%
  #.[, `:=`(mean_i = NULL, mean_j = NULL)] %>%
  .[]

# 5. resampled data ============================================================
fpm_resampled_summary <- fpm_resampled[, list(
  mean = mean(value),
  SD = sd(value),
  var = var(value),
  #CV = sd(value) / mean(value),
  Q_.025 = quantile(value, 0.025),
  Q_.25 = quantile(value, 0.25),
  Q_.5 = quantile(value, 0.5),
  Q_.75 = quantile(value, 0.75),
  Q_.975 = quantile(value, 0.975)
), by = .(variable, case, j)]

fpn_resampled_summary <- fpn_resampled[, list(
  mean = mean(value, na.rm = TRUE),
  SD = sd(value),
  var = var(value),
  #CV = sd(value) / mean(value),
  Q_.025 = quantile(value, 0.025),
  Q_.25 = quantile(value, 0.25),
  Q_.5 = quantile(value, 0.5),
  Q_.75 = quantile(value, 0.75),
  Q_.975 = quantile(value, 0.975)
), by = .(variable, case)]

fpm_resampled_summary[variable == 'Fmat', variable := 'orig']
fpm_resampled_summary[variable != 'Fmat', variable := gsub('Fmat_resampled_', '',
                                                           variable)]

fpn_resampled_summary[variable == 'Fmat', variable := 'orig']
fpn_resampled_summary[variable != 'Fmat', variable := gsub('Fmat_resampled_', '',
                                                           variable)]


# fpm_resampled_summary <- rbindlist(list(fpm_summary[, .(case, j, mean, sd, cv, CI25, CI975, median)],
#                                         fpm_resampled_summary), fill = TRUE)
# fpm_resampled_summary[is.na(type), type := 'orig']
fpm_resampled_summary <- melt(fpm_resampled_summary,
                              id.vars = c('case', 'variable', 'j'),
                              variable.name = 'measure')
fpm_resampled_summary <- merge(fpm_resampled_summary[variable != 'orig'],
                               fpm_resampled_summary[variable == 'orig'],
                               by = c('case', 'j', 'measure'),
                               suffixes = c('', '_orig'), sort = FALSE)

fpm_resampled_summary <- merge(fpm_resampled_summary,
                               fpm_summary[, .(j, case, abs_mean, abs_mean_rel,
                                               CodeNr)],
                               by = c('j', 'case'), sort = FALSE)

fpm_resampled_summary[, absdif := value - value_orig]
fpm_resampled_summary[, reldif := absdif / value_orig]
fpm_resampled_summary[, perc_dev := (value / value_orig) - 1]
fpm_resampled_summary[, ratio := (value / value_orig)]

fpm_resampled_summary[, logratio := log10(value /value_orig)]


# fpm_resampled_summary <- fpm_resampled_summary %>%
#   melt(id.vars = c('case', 'variable', 'j'), variable.name = 'measure') %>%
#   dcast(case + j + measure ~ variable, value.var = 'value')
# fpm_resampled_summary[, absdif_mv := mv - orig]
# fpm_resampled_summary[, reldif_mv := absdif_mv / orig]
#
# fpm_resampled_summary[, absdif_uv := uv - orig]
# fpm_resampled_summary[, reldif_uv := absdif_uv / orig]


#fpm_resampled_summary[, dev_mv := (mv / orig) - 1]
#fpm_resampled_summary[, dev_uv := (uv / orig) - 1]





# fpn_resampled_summary <- rbindlist(list(fpn_summary[, .(case, mean, sd, cv, CI25, CI975, median)],
#                                         fpn_resampled_summary), fill = TRUE)
# fpn_resampled_summary[is.na(type), type := 'orig']
fpn_resampled_summary <- fpn_resampled_summary %>%
  melt(id.vars = c('case', 'variable'), variable.name = 'measure') %>%
  dcast(case + measure ~ variable, value.var = 'value')

# fpn_resampled_sample <- rbindlist(list(fpn_sample, fpn_resampled),
#                                   fill = TRUE)
# fpn_resampled_sample[is.na(type), type := 'orig']

fpn_resampled[variable == 'Fmat', label := 'original']
fpn_resampled[variable != 'Fmat', label := str_replace(variable, 'Fmat_resampled_', '')]


############################################################################## #
##### Plots #############################################################
############################################################################## #


# 1. National footprints =======================================================
# _a) histogram =============
ggplot(fpn_sample, aes(x = value, fill = as.factor(case),
                       col = as.factor(case))) +
  geom_histogram(position = 'identity', alpha = 0.1, linewidth = 1.3) +
  #  facet_wrap(~case) +
  scale_fill_colorblind7() +
  scale_color_colorblind7()




# _b) ridge plot ============================================
# this one!
linewidth = 0.4
(p1 <- ggplot(fpn_sample, aes(x = value, y = as.factor(case),
                              fill = as.factor(case),
                              col = as.factor(case) )) +
    geom_vline(data = fpn_summary,
               aes(xintercept = median, col = as.factor(case)),
               linetype = 'solid', linewidth = linewidth) +
    geom_vline(data = fpn_summary,
               aes(xintercept = CI25, col = as.factor(case)),
               linetype = 'dashed', linewidth = linewidth) +
    geom_vline(data = fpn_summary,
               aes(xintercept = CI975, col = as.factor(case)),
               linetype = 'dashed', linewidth = linewidth) +
    # stat_density_ridges(alpha=0, quantile_lines = TRUE,
    #                     quantiles = c(0.025, 0.5, 0.975), linetype = 'dashed') +
    # geom_density_ridges(alpha=0, linewidth = 1.2, col = 'white') +
    geom_density_ridges(alpha=0.6, stat="binline", bins = 30, linewidth = 1.2) +
    scale_fill_colorblind() +
    scale_color_colorblind() +
    ylab('Case') +
    xlab('Emissions [Gg]') +
    theme_ridges(font_size = 11) +
    theme(legend.position = 'none')
)

ggsave(filename = file.path(path2plots, "fp_national_ridge.pdf"),
       height = 5, width = 7)

# 2. Multiplier =======================================================

# _a) ridge plot ==============================================
linewidth = 0.4
ggplot(fpm_summary, aes(x = cv, y = as.factor(case),
                        fill = as.factor(case),
                        col = as.factor(case) )) +
  geom_vline(data = fpm_cv_summary,
             aes(xintercept = median, col = as.factor(case)),
             linetype = 'solid', linewidth = linewidth) +
  geom_vline(data = fpm_cv_summary,
             aes(xintercept = CI25, col = as.factor(case)),
             linetype = 'dashed', linewidth = linewidth) +
  geom_vline(data = fpm_cv_summary,
             aes(xintercept = CI975, col = as.factor(case)),
             linetype = 'dashed', linewidth = linewidth) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 40, linewidth = 1.2) +
  scale_x_log10() +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  theme_ridges(font_size = 11) +
  theme(legend.position = 'none')

# _b) violin plot ==========================================================
# this one!
(p2 <- ggplot(fpm_summary, aes(y = cv, x = as.factor(case))) +
   #scale_y_log10() +
   stat_summary(fun.data =  median.quartile, geom = 'crossbar', width = .2,
                fill = 'white') +
   geom_jitter(aes(size = mean_rel), shape=16, alpha = 0.4, col = 'grey50',
               position=position_jitter(0.1, seed = 1)) +
   stat_summary(fun.data =  percentile, geom = 'errorbar', width = .2,
                size = 1, col = my_cols[3]) +
   geom_violin(trim = TRUE, alpha = 0.2, fill = my_cols[7], col = NA) +
   stat_summary(fun.data =  median.quartile,
                geom = 'crossbar', width = .2, size = 0.7,
                alpha=0.4, fill = my_cols[3], col = my_cols[3]) +
   # facet_zoom(ylim = c(0, 0.5), zoom.size = 1) +
   ylab('CV = sd / mean')  +
   xlab('Case') +
   labs(size = 'Share of total', alpha = 'Share of total') +
   theme(legend.position = 'bottom',
         strip.background  = element_rect(fill = 'grey90', colour = 'grey70'))
)

ggsave(filename = file.path(path2plots, "fp_multiplier_violin.pdf"),
       height = 6, width = 10)


wrap_plots(p2, p1, widths = c(2,1)) +
  plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(path2plots, "figure5.pdf"),
       height = 6, width = 10)

# _c) Paper 2 type plot  ============================================================
setorder(fpm_summary, case, mean)
fpm_summary[, right := cumsum(mean)/ sum(mean), by = .(case)]
fpm_summary[, left := right - (mean / sum(mean)), by = .(case)]


ggplot(fpm_summary, aes(xmin = left, xmax = right, ymax = (CI975 / mean) - 1,
                        ymin = (CI25 / mean) - 1, y = (median / mean) - 1)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_rect(col = 'grey50', size = 0.2, fill = my_cols[2], alpha = 0.7) +
  geom_text_repel(aes(label = j, y = (CI975 / mean) - 1, x = left + (right - left) / 2),
                  size = 2.3, position = position_nudge_repel(y = 0.1),
                  col = my_cols[4], max.overlaps = 15) +
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)+
  facet_wrap(~case, scales = 'free') +
  ylab('% deviation from sample mean') +
  xlab('Share of total emissions') +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(legend.position = 'bottom',
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background  = element_rect(fill = 'grey90', colour = NA))

ggplot(fpm_summary, aes(xmin = left, xmax = right, ymax = cv,
                        ymin = 0, y = (median / mean) - 1)) +
  #geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_rect(col = 'grey50', size = 0.2, fill = my_cols[2], alpha = 0.7) +
  scale_shape_manual(name = "Legend", values = c('Median' = 3, 'CV' = 4, 'EXIOBASE' = 17)) +
  scale_color_manual(name = "Legend", values = c('Median' = 1, 'CV' = my_cols[7], 'EXIOBASE' = my_cols[3])) +
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)+
  facet_wrap(~case) +
  ylab('% deviation from sample mean') +
  xlab('Share of total emissions') +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(legend.position = 'bottom',
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background  = element_rect(fill = 'grey90', colour = NA))


# _d) full sample ===============================================

ggplot(fpm_sample[case == 1], aes(x = value, y = as.factor(j))) +

  # stat_density_ridges(alpha=0, quantile_lines = TRUE,
  #                     quantiles = c(0.025, 0.5, 0.975), linetype = 'dashed') +
  # geom_density_ridges(alpha=0, linewidth = 1.2, col = 'white') +
  geom_density_ridges(alpha=0.6,linewidth = 0.2) +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  scale_x_log10() +
  facet_wrap(~case, scales = 'free') +
  theme_ridges() +
  theme(legend.position = 'none')

# 3. AEA =============================================================

# _a) ridge plot ============================================
linewidth = 0.4
ggplot(aea_summary[], aes(x = cv, y = as.factor(case),
                          fill = as.factor(case),
                          col = as.factor(case) )) +
  # geom_vline(data = fpm_cv_summary,
  #            aes(xintercept = median, col = as.factor(case)),
  #            linetype = 'solid', linewidth = linewidth) +
  # geom_vline(data = fpm_cv_summary,
  #            aes(xintercept = CI25, col = as.factor(case)),
  #            linetype = 'dashed', linewidth = linewidth) +
  # geom_vline(data = fpm_cv_summary,
  #            aes(xintercept = CI975, col = as.factor(case)),
  #            linetype = 'dashed', linewidth = linewidth) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 40, linewidth = 1.2) +
  scale_x_log10() +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  theme_ridges() +
  theme(legend.position = 'none')

# _b) point plot =========================================
linewidth = 0.4
ggplot(aea_summary[], aes(x = cv, y = mean,
                          fill = as.factor(case),
                          col = as.factor(case) )) +
  geom_point(alpha = 0.4) +
  scale_y_log10() +
  facet_zoom(xlim = c(0, 2), ylim = c(0, aea_summary[, log10(max(mean))]), zoom.size = 2) +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  theme(legend.position = 'none')



# _c) Paper 2 type plot  ============================================================
setorder(aea_summary, case, mean)
aea_summary[, right := cumsum(mean)/ sum(mean), by = .(case)]
aea_summary[, left := right - (mean / sum(mean)), by = .(case)]


ggplot(aea_summary, aes(xmin = left, xmax = right, ymax = (CI97.5 / mean) - 1,
                        ymin = (CI2.5 / mean) - 1, y = (median / mean) - 1)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_rect(col = 'grey50', size = 0.2, fill = my_cols[2], alpha = 0.7) +
  geom_text_repel(aes(label = y, y = (CI97.5 / mean) - 1, x = left + (right - left) / 2),
                  size = 2.3, position = position_nudge_repel(y = 0.1),
                  col = my_cols[4], max.overlaps = 15) +
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)+
  facet_wrap(~case, scales = 'fixed') +
  ylab('% deviation from sample mean') +
  xlab('Share of total emissions') +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(legend.position = 'bottom',
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background  = element_rect(fill = 'grey90', colour = NA))

ggplot(aea_summary[left > 0.05], aes(xmin = left, xmax = right, ymax = cv,
                                     ymin = 0)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_rect(col = 'grey50', size = 0.2, fill = my_cols[2], alpha = 0.7) +
  geom_text_repel(aes(label = y, y = cv, x = left + (right - left) / 2),
                  size = 2.3, position = position_nudge_repel(y = 0.1),
                  col = my_cols[4], max.overlaps = 15) +
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)+
  facet_wrap(~case, scales = 'fixed', ncol = 2) +
  ylab('CV = SD / mean') +
  xlab('Share of total emissions') +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(legend.position = 'bottom',
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background  = element_rect(fill = 'grey90', colour = NA))





# 3. Correlations =============================================================

# _a) footprints ==============================================================
# __i. violin box plot ===========================
# this one!
ggplot(fpm_cor_sample, aes(y = value, x = as.factor(case))) +
  #scale_y_log10() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(fun.data =  median.quartile, geom = 'crossbar', width = .2,
               fill = 'white') +
  # geom_jitter(shape=16, alpha = 0.2, col = 'grey50', size = 0.6,
  #             position=position_jitter(0.1, seed = 1)) +
  stat_summary(fun.data =  percentile, geom = 'errorbar', width = .2,
               size = 1, col = my_cols[3]) +
  geom_violin(trim = TRUE, alpha = 0.2, fill = my_cols[7], col = NA) +
  stat_summary(fun.data =  median.quartile,
               geom = 'crossbar', width = .2, size = 0.7,
               alpha=0.4, fill = my_cols[3], col = my_cols[3]) +
  #facet_zoom(ylim = c(0, 1)) +
  ylab('rho')  +
  xlab('Case') +
  labs(size = 'Share of total', alpha = 'Share of total') +
  theme(legend.position = 'bottom')

ggsave(filename = file.path(path2plots, "fp_multiplier_cor.pdf"),
       height = 5, width = 7)


# __ii. ridge plot ========================================
linewidth = 0.4
ggplot(fpm_cor_sample, aes(x = value, y = as.factor(case),
                           fill = as.factor(case),
                           col = as.factor(case) )) +
  # geom_vline(data = fpm_cor_summary,
  #            aes(xintercept = median, col = as.factor(case)),
  #            linetype = 'solid', linewidth = linewidth) +
  # geom_vline(data = fpm_cor_summary,
  #            aes(xintercept = CI25, col = as.factor(case)),
  #            linetype = 'dashed', linewidth = linewidth) +
  # geom_vline(data = fpm_cor_summary,
  #            aes(xintercept = CI975, col = as.factor(case)),
  #            linetype = 'dashed', linewidth = linewidth) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 40, linewidth = 1.2) +
  #scale_x_log10() +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  theme_ridges() +
  theme(legend.position = 'none')


# __ii. xy plot ==============================================
(p1 <- ggplot(fpm_cor_sample[abs(value) > 0.1 & abs_mean_i > 1E-2
                             & abs_mean_j > 1E-2],
              aes(x = abs_mean_i, y = abs_mean_j, col = value)) +
   geom_point(shape= 16) +
   geom_point(data = fpm_cor_sample[abs(value) > 0.1 & abs_mean_i > 1E-2
                                    & abs_mean_j > 1E-2 & value < 0]) +
   # scale_x_log10() +
   # scale_y_log10() +
   xlab(bquote(~mean[i])) +
   ylab(bquote(~mean[j])) +
   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))) +
   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))) +
   annotation_logticks(color = 'grey60', alpha = 0.4) +
   facet_wrap(~case, scales = 'fixed', nrow = 1) +
   theme(aspect.ratio=1, legend.position = 'bottom') +
   scale_color_gradient2(limits = c(-1, 1))
)

ggsave(filename = file.path(path2plots, "fpm_cor_xy.pdf"),
       height = 5, width = 7)

(ggplot(fpm_cor_sample[(value) > 0.1 ],
        aes(x = cv_i, y = cv_j, col = value, alpha = abs(value),
            size = abs(value))) +
    geom_point() +
    geom_point(data = fpm_cor_sample[value < -0.1 ], inherit.aes = TRUE) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~case, scales = 'fixed', nrow = 1) +
    theme(aspect.ratio=1, legend.position = 'bottom') +
    scale_color_gradient2()
)
# _b) AEA ===================================================
library(scales)
#xy plot
(p2 <- ggplot(aea_cor2[abs(value) > 0.1],
              aes(x = mean_i, y = mean_j,col = value)) +
    geom_point(shape = 16) +
    geom_point(data = aea_cor2[(value) < -0.1], shape = 16) +
    #scale_x_log10() +
    #scale_y_log10() +
    # labs(x = bquote(~a_i)) +
    xlab(bquote(~mean[i])) +
    ylab(bquote(~mean[j])) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(color = 'grey60', alpha = 0.4) +
    facet_wrap(~case, scales = 'fixed', nrow = 1) +
    theme(aspect.ratio=1, legend.position = 'none') +
    scale_color_gradient2(limits = c(-1, 1))
)

ggsave(filename = file.path(path2plots, "aea_cor_xy.pdf"),
       height = 5, width = 7)

wrap_plots(p2, p1, nrow = 2, guides = 'auto') +
  plot_annotation(tag_levels = 'A')

ggsave(filename = file.path(path2plots, "cor_xy_combined.pdf"),
       height = 5, width = 7)




linewidth = 0.4
ggplot(aea_cor2, aes(x = value, y = as.factor(case),
                     fill = as.factor(case),
                     col = as.factor(case) )) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 80, linewidth = 1.2) +
  #scale_x_log10() +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  theme_ridges() +
  theme(legend.position = 'none')

ggplot(aea_cor2[case > 2], aes(x = mean, y = value, col = as.factor(case))) +
  geom_point(alpha = 0.3) +
  scale_x_log10() +
  scale_color_colorblind()

setorder(aea_cor2, case, mean)
aea_cor2[, right := cumsum(mean)/ sum(mean), by = .(case)]
aea_cor2[, left := right - (mean / sum(mean)), by = .(case)]

ggplot(aea_cor2, aes(xmin = left, xmax = right, ymax = value,
                     ymin = 0)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_rect(col = 'grey50', size = 0.2, fill = my_cols[2], alpha = 0.7) +
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(labels = scales::percent)+
  facet_wrap(~case, scales = 'fixed', ncol = 2) +
  ylab('rho') +
  xlab('Share of total emissions') +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(legend.position = 'bottom',
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background  = element_rect(fill = 'grey90', colour = NA))


# violin box plot
ggplot(aea_cor2, aes(y = value, x = as.factor(case))) +
  #scale_y_log10() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_summary(fun.data =  median.quartile, geom = 'crossbar', width = .2,
               fill = 'white') +
  geom_jitter(shape=16, alpha = 0.2, col = 'grey50', size = 0.6,
              position=position_jitter(0.1, seed = 1)) +
  stat_summary(fun.data =  percentile, geom = 'errorbar', width = .2,
               size = 1, col = my_cols[3]) +
  geom_violin(trim = TRUE, alpha = 0.2, fill = my_cols[7], col = NA) +
  stat_summary(fun.data =  median.quartile,
               geom = 'crossbar', width = .2, size = 0.7,
               alpha=0.4, fill = my_cols[3], col = my_cols[3]) +
  #facet_zoom(ylim = c(0, 1)) +
  ylab('rho')  +
  xlab('Case') +
  labs(size = 'Share of total', alpha = 'Share of total') +
  theme(legend.position = 'bottom')



# 4. Resample =================================================================

# _a) Multiplier ==================================================
# __i. boxplots ===================


# __ii. Only case 3 ==================================================

plot_data <- fpm_resampled_summary %>%
  .[measure %in% c('var')] %>% #, 'Q_.25', 'Q_.75'
  .[variable %in% c('mv_gamma', 'uv_gamma')] %>%
  .[case == 3]

mylims <- range(with(plot_data, c(value_orig, value)), na.rm=TRUE)

plot_data[(logratio == max(logratio)
           | logratio == min(logratio)),
          label := CodeNr]
plot_data[!is.na(label)]



(p_var <- plot_data %>%
    .[abs_mean_rel > 1E-8] %>%
    ggplot(aes(x = abs_mean_rel,
               y = log2(value / value_orig),#reldif,
               col = variable)) +
    geom_abline(slope = 0, intercept = 0, linetype = 'dashed') +
    geom_point(alpha = 0.6, shape = 16, show.legend = FALSE) +
    geom_point(data = plot_data[!is.na(label)], alpha = 1, size = 1.5,
               show.legend = FALSE) +
    geom_text_repel(aes(label = label), col = 'black') +
    scale_color_colorblind7() +
    scale_y_continuous(labels = function(breaks) (round(2^(breaks), 1)),
                       breaks = scales::breaks_extended(n = 9)) +

    facet_wrap(~measure, scales = 'free_y') +
    #  tune::coord_obs_pred(xlim = mylims, ylim = mylims) +
    ylab('gamma_sample / original_sample')  +
    xlab('Share of total impact') +
    scale_x_log10() +
    #  scale_y_log10() +
    labs(col = NULL, fill = NULL) +
    theme(axis.title.y = element_blank())
)

(p_mean <- fpm_resampled_summary %>%
    .[measure %in% c('mean')] %>% #, 'Q_.25', 'Q_.75'
    .[variable %in% c('mv_gamma', 'uv_gamma')] %>%
    .[case == 3] %>%
    .[abs_mean_rel > 1E-8] %>%
    ggplot(aes(x = abs_mean_rel,
               y = log2(value / value_orig),#reldif,
               col = variable)) +
    geom_abline(slope = 0, intercept = 0, linetype = 'dashed') +
    geom_point(alpha = 0.6, shape = 16, show.legend = FALSE) +
    scale_color_colorblind7() +
    scale_y_continuous(labels = function(breaks) (round(2^(breaks), 3)),
                       breaks = scales::breaks_extended(n = 9)) +

    facet_wrap(~measure, scales = 'free_y') +
    #  tune::coord_obs_pred(xlim = mylims, ylim = mylims) +
    ylab('gamma_sample / original_sample')  +
    xlab('Share of total impact') +
    scale_x_log10() +
    #  scale_y_log10() +
    labs(col = NULL, fill = NULL) +
    theme(legend.position = 'none')
)


ggsave(filename = file.path(path2plots, "fp_multiplier_resampled.pdf"),
       height = 5, width = 7)




# __iii. compare gamma vs lnorm ==================================
fpm_resampled_summary$variable %>% unique


plot_data <- fpm_resampled_summary %>%
  .[measure %in% c('var', 'mean')] %>% #, 'Q_.25', 'Q_.75'
  # .[variable %in% c('mv_gamma', 'uv_gamma', 'mv', 'uv')] %>%
  .[case == 3]

plot_data[grepl('^mv', variable), type := 'multivariate']
plot_data[grepl('^uv', variable), type := 'univariate']

plot_data[variable %in% c('mv', 'uv'), dist_type := 'Lognormal']
plot_data[variable %in% c('mv_gamma', 'uv_gamma'), dist_type := 'Gamma']
plot_data[grepl('truncnorm', variable), dist_type := 'Truncated\nNormal']

plist <- list()
for (imeasure in c('var', 'mean')) {
  plist[[imeasure]] <- list()
  for (itype in c('univariate', 'multivariate')) {
    mylims <- range(with(plot_data, c(value_orig, value)), na.rm=TRUE)


    (plist[[imeasure]][[itype]] <- plot_data[measure == imeasure & type == itype] %>%
        .[abs_mean_rel > 1E-8] %>%
        ggplot(aes(x = dist_type,
                   y = log2(value / value_orig),#reldif,
                   col = dist_type)) +
        geom_abline(slope = 0, intercept = 0, linetype = 'dashed') +
        #geom_point(alpha = 0.6, shape = 16, show.legend = TRUE) +
        geom_jitter(alpha = 0.3, shape = 16) +
        geom_boxplot(outlier.alpha = 0, fill = NA) +
        scale_color_colorblind7() +
        scale_y_continuous(
          labels = function(breaks) (round(2^(breaks),
                                           ifelse(imeasure == 'var', 2, 3))),
          breaks = scales::breaks_extended(n = 9)) +

        #facet_wrap(~measure, scales = 'free_y') +
        #  tune::coord_obs_pred(xlim = mylims, ylim = mylims) +
        ylab('ratio: resampled / original')  +
        xlab('Type of distribution') +
        #     scale_x_log10() +
        ggtitle(paste0(imeasure, ' ', itype)) +
        theme(legend.position = 'none')
    )
  }
}


wrap_plots(
  plist$var$multivariate + ggtitle("Multivariate sampling"),
  plist$var$univariate + ggtitle("Univariate sampling"),
  guides = 'collect'
) & theme(legend.position = 'none') &
  labs(col = 'Type of distribution')

plist$var$multivariate +
  ggtitle("") +
  ylab('ratio of the variances: var(resampled) / var(original)')

ggsave(filename = file.path(path2plots, "figureS2.pdf"),
       height = 5, width = 5)


# __iv. all cases ==================================










# _b) National fps ==================================================

# __i. Gamma distribtuion ===========================================
fpn_resampled$label %>% unique
fpn_resampled[, label := factor(label, levels = c('original',
                                                  'mv_gamma',
                                                  'uv_gamma'))]

# this one!
linewidth = 0.5

(p3 <- ggplot(fpn_resampled[label %in% c('original',
                                         'mv_gamma', 'uv_gamma') &
                              case == 3],
              aes(x = value, #y = as.factor(case),
                  fill = label,
                  col = label )) +
    geom_histogram(alpha=0.3, position = 'identity', linewidth = linewidth) +
    scale_fill_colorblind(labels = c("Original\nsample",
                                     "Multivariate\nsample",
                                     "Univariate\nsample")) +
    scale_color_colorblind(labels = c("Original\nsample",
                                      "Multivariate\nsample",
                                      "Univariate\nsample")) +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    #scale_fill_discrete(name = "", ) +
    #ylab('Case') +
    xlab('Emissions [Gg]') +
    labs(col = NULL, fill = NULL) +
    #facet_wrap(~case, labeller = 'label_both') +
    # theme_ridges(font_size = 11) +
    theme(legend.position = 'bottom')
)


ggsave(p3, filename = file.path(path2plots, "fp_national_resampled.pdf"),
       height = 4, width = 7)

wrap_plots(wrap_plots(p_mean, p_var), p3, guides = 'collect', widths = c(2,1)) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'bottom')

ggsave(filename = file.path(path2plots, "figure6.pdf"),
       height = 4, width = 9)

# __ii. Compare gamma vs. lnorm ======================================

plot_data <- fpn_resampled[label %in% c('original',
                                        'mv_gamma', 'uv_gamma',
                                        'mv', 'uv') &
                             case == 3]

plot_data$label %>% unique

plot_data[grepl('^mv', label), type := 'multivariate']
plot_data[grepl('^uv', label), type := 'univariate']
plot_data[label == 'original', type := label]

plot_data[label %in% c('mv', 'uv'), dist_type := 'Lognormal']
plot_data[label %in% c('mv_gamma', 'uv_gamma'), dist_type := 'Gamma']
plot_data[label == 'original', dist_type := label]



plot_data[, dist_type := factor(dist_type, levels = c('original',
                                                      'Gamma',
                                                      'Lognormal'))]

linewidth = 0.5

plist2 <- list()
for (itype in c('univariate', 'multivariate')) {
  (plist2[[itype]] <- ggplot(plot_data[type %in% c(itype, 'original')],
                             aes(x = value, #y = as.factor(case),
                                 fill = dist_type,
                                 col = dist_type )) +
     geom_histogram(alpha=0.3, position = 'identity', linewidth = linewidth) +
     scale_fill_colorblind(labels = c("Original\nsample",
                                      "Gamma",
                                      "Lognormal")) +
     scale_color_colorblind(labels = c("Original\nsample",
                                       "Gamma",
                                       "Lognormal")) +
     scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
     #scale_fill_discrete(name = "", ) +
     #ylab('Case') +
     xlab('Emissions [Gg]') +
     labs(col = NULL, fill = NULL) +
     #facet_wrap(~case, labeller = 'label_both') +
     # theme_ridges(font_size = 11) +
     theme(legend.position = 'bottom')
  )


}
wrap_plots(plist2)

ggsave(p3, filename = file.path(path2plots, "fp_national_resampled.pdf"),
       height = 4, width = 7)

wrap_plots(wrap_plots(p_mean, p_var), p3, guides = 'collect', widths = c(2,1)) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'bottom')

ggsave(filename = file.path(path2plots, "fp_combined_resampled.pdf"),
       height = 4, width = 9)


# ____archive =========================

ggplot(fpn_resampled_sample[case == 3],
       aes(x = value,
           fill = type,
           col = type )) +
  # geom_vline(data = fpn_summary,
  #            aes(xintercept = median, col = as.factor(case)),
  #            linetype = 'solid', linewidth = linewidth) +
  # geom_vline(data = fpn_summary,
  #            aes(xintercept = CI25, col = as.factor(case)),
  #            linetype = 'dashed', linewidth = linewidth) +
  # geom_vline(data = fpn_summary,
  #            aes(xintercept = CI975, col = as.factor(case)),
  #            linetype = 'dashed', linewidth = linewidth) +
  # # stat_density_ridges(alpha=0, quantile_lines = TRUE,
  #                     quantiles = c(0.025, 0.5, 0.975), linetype = 'dashed') +
  # geom_density_ridges(alpha=0, linewidth = 1.2, col = 'white') +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  scale_x_log10() +
  ylab('Case') +
  xlab('Emissions [Gg]') +
  theme_ridges() +
  theme(legend.position = 'bottom')





# _c) Identify multipliers with largest dicrepencies when resampling============


fpm_resampled_summary %>%
  .[measure %in% c('SD')]%>%
  .[variable %in% c('uv_gamma')] %>%
  .[case == 3 & ratio > 1.5]

fpm_resampled_summary %>%
  .[measure %in% c('SD')]%>%
  .[variable %in% c('uv_gamma')] %>%
  .[ratio == min(ratio)]

fpm_resampled_summary %>%
  .[measure %in% c('var')]%>%
  .[variable %in% c('uv_gamma')] %>%
  .[ratio > -1E-0 & ratio < 1E-0]


fpm_resampled_summary %>%
  .[measure %in% c('SD')]%>%
  .[variable %in% c('uv_gamma')] %>%
  .[case == 2] %>%
  ggplot(aes(x = value_orig, y = logratio)) +
  # scale_x_log10() +
  geom_point()
## mmh TODO: all sectors with a large logratio (overest with UV gamma) are very small. might be just an artefact??


# sectors with high ratio
# Transmission of electricity     DE   923
# Distribution and trade of electricity     DE   924
# Pulp     DE   867

# sector with low ratio
# Manufacture of wearing apparel; dressing and dyeing of fur (18)     DE   863
indices[col == 863]

i_overest <- 923
i_underest <- 863
# look at supply chain of those sectors
Lmat <- readRDS(file.path(path2output, 'prepare_EXIOBASE_L.RData'))
xvec <- readRDS(file.path(path2output, 'prepare_EXIOBASE_x.RData'))

irows <- indices[region == 'DE']$col

plot_cormat_for_sector_inputs <- function(i) {
  L_dt <- data.table(
    id = irows,
    industry_code = indices[region == 'DE']$CodeNr,
    L = Lmat[irows,i],
    x = xvec[irows]
  )


  dt <- merge(L_dt, aea_raw[case == 3]$sample[[1]],
              by.x = 'industry_code', by.y = 'y', all = TRUE)
  dt <- dt[industry_code != 'y01']

  # calculate S = F / x
  dt[, sample_S := pmap(list(Fmat = sample, x = x),
                        function(Fmat, x) Fmat / x)]

  # calculate multipliers M = S*L
  dt[, sample_M := pmap(list(Smat = sample_S, L = L),
                        function(Smat, L) Smat * L)]

  # throw out zero inputs
  dt <- dt[sapply(sample_M, length) > 0]

  # check correlations among inputs
  inputs_matrix <- do.call('cbind', dt$sample_M)
  colnames(inputs_matrix) <- dt$industry_code

  # test if everything worked
  #rowSums(inputs_matrix) %>% mean
  #fpm_summary[j == i_underest]

  cor_dt <- cordt(inputs_matrix)

  mean_dt <- data.table(mean = colMeans(inputs_matrix),
                        id = colnames(inputs_matrix))

  cor_dt2 <- merge(cor_dt, mean_dt, by.x = 'i', by.y = 'id')
  cor_dt3 <- merge(cor_dt2, mean_dt, by.x = 'j', by.y = 'id', suffixes = c('_x', '_y'))

  (ggplot(cor_dt3[mean_x>1E-8 & mean_y > 1E-8 &
                    abs(value) > 0.3],
          aes(x = mean_x, y = mean_y, col = value)) +
      geom_point(shape= 16, size = 2.4) +
      # scale_x_log10() +
      # scale_y_log10() +
      xlab(bquote(~mean[i])) +
      ylab(bquote(~mean[j])) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      annotation_logticks(color = 'grey60', alpha = 0.4) +
      #facet_wrap(~case, scales = 'fixed', nrow = 1) +
      theme(aspect.ratio=1, legend.position = 'bottom') +
      labs(color = 'Correlation coefficient') +
      ggtitle(paste(unlist(indices[col == i, .(CodeNr, sector)]),
                    collapse = ': ')) +
      scale_color_gradient2(limits = c(-1, 1))
  )

}
p1 <- plot_cormat_for_sector_inputs(i_overest)
p2 <- plot_cormat_for_sector_inputs(i_underest)

plot_cormat_for_sector_inputs(817)


wrap_plots(p1, p2, guides = 'collect') &
  theme(legend.position = 'bottom')

ggsave(filename = file.path(path2plots, "figureS3.pdf"),
       height = 5, width = 7)


cor_dt3[mean_y > 1E-3 & abs(value) > 0.3 & value > 0] %>% setorder(mean_x) %>% .[]
sum(mean_dt$mean)
# 90% of inputs come from Electr. Coal, which is negatively correlated to



############################################################################## #
##### save results #############################################################
############################################################################## #


# THE END ---------------------------------------------------------------------
