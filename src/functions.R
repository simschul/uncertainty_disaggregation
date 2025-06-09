library(patchwork)
#library(RinR)
library(grid)
library(ggthemes)

compare_samples <- function(sample_orig, sample_new, ks_test = FALSE, print = FALSE,
                            type = 'hist', log = FALSE, return = 'both') {

  if (isTRUE(print)) {
    print('Summary sample_orig:')
    print(summary(sample_orig))

    print('Summary sample_new:')
    print(summary(sample_new))
  }

  sample_orig_agg <- rowSums(sample_orig)
  sample_new_agg <- rowSums(sample_new)

  sample_orig_agg <- reshape2::melt(sample_orig_agg) %>%
    suppressWarnings %>%
    as.data.table(keep.rownames = 'Var1') %>%
    .[, Var1 := as.integer(Var1)]
  sample_new_agg <- reshape2::melt(sample_new_agg) %>%
    suppressWarnings %>%
    as.data.table(keep.rownames = 'Var1') %>%
    .[, Var1 := as.integer(Var1)]

  sample_orig <- sample_orig %>%
    melt %>%
    suppressWarnings %>%
    as.data.table

  sample_new <- sample_new %>%
    melt %>%
    suppressWarnings %>%
    as.data.table

  if (type %in% c('hist', 'ks')) {

    data <- merge(sample_orig, sample_new, by = c('Var1', 'Var2'),
                  suffixes = c("_orig", '_new'))

    if (isTRUE(ks_test)) {
      p_values <- data[, list(p_value = paste0('p = ',
                                               round(ks.test(value_new, value_orig)$p.value, 3)),
                              x = max(c(value_new, value_orig)) - (0.1 * (max(c(value_new, value_orig)) - min(c(value_new, value_orig)))),
                              y = nrow(sample_orig)/100,
                              variable = 'value_orig'),
                       by = .(Var2)]

    }



    data2 <- melt(data, id.vars =  c('Var1', 'Var2'),
                  measure.vars = c('value_orig', 'value_new'))
    # data_agg <- data2[, list(value = sum(value)), by = .(variable, Var1)]
    data_agg <- rbindlist(list('value_orig' = sample_orig_agg,
                               'value_new' = sample_new_agg), idcol = 'variable')
    #data_combined <- merge(data2, data_agg, by = c('Var1', 'variable'),
    #              )

    data_agg[, Var2 := as.integer(0)]
    data_combined <- rbindlist(list(disagg = data2, agg = data_agg),
                               idcol = 'type', use.names = TRUE)

    # plot
    if (type == 'ks') {
      # ggplot(data_combined, aes(x = value, linetype = variable,
      #                   col = as.factor(Var2))) +
      #   stat_ecdf(geom = "step", pad = FALSE, linewidth = 1, alpha = 0.7) +
      #   #facet_wrap(~Var2, scales = 'free_x') +
      #   scale_color_colorblind() +
      #   theme(legend.position = 'right')


      p1 <- ggplot(data2, aes(x = value, linetype = variable,
                              col = as.factor(Var2))) +
        stat_ecdf(geom = "step", pad = FALSE, linewidth = 1, alpha = 0.7) +
        #facet_wrap(~Var2, scales = 'free_x') +
        scale_color_colorblind7()


      p2 <- ggplot(data_agg, aes(x = value, linetype = variable)) +
        stat_ecdf(geom = "step", pad = FALSE, linewidth = 1, alpha = 0.7) +
        scale_color_colorblind7() +
        theme(axis.title.x=element_blank())
      #p1 / p2
      #cowplot::plot_grid(p2, p1, ncol = 1)
      p_combined <- p2 / p1 & theme(legend.position = "bottom")
      p_combined + plot_layout(guides = "collect")

    } else if (type == 'hist') {
      data2[variable == 'value_orig', colour := as.factor(1)]
      data2[variable != 'value_orig', colour := as.factor(Var2+1)]

      data_agg[variable == 'value_orig', colour := as.factor(2)]
      data_agg[variable != 'value_orig', colour := as.factor(1)]

      p1 <- ggplot(data2, aes(x = value, col = colour, fill = colour)) +
        geom_histogram(alpha = 0.3, position = 'identity') +
        facet_wrap(~Var2, scales = 'free_x') +
        #geom_text(data = p_values, aes(x=x, y=y, label = p_value)) +
        scale_fill_colorblind3() +
        scale_color_colorblind3() +
        theme(legend.position = 'none',
              strip.background = element_blank(),
              strip.text.x = element_blank())

      if (isTRUE(log)) {
        p1 <- p1 + scale_x_log10()
      }

      if (isTRUE(ks_test)) {
        p1 <- p1 + geom_text(data = p_values, aes(x=x, y=y, label = p_value))
        p_value_agg <- data.table(
          variable = 'value_orig',
          p_value = paste0('p = ',
                           round(ks.test(data_agg[variable == 'value_new']$value,
                                         data_agg[variable == 'value_orig']$value)$p.value, 3)),
          x = max(c(data_agg$value)) - (0.1 * (max(c(data_agg$value)) - min(c(data_agg$value)))),
          y = nrow(sample_orig)/100
        )

      }

      p2 <- ggplot(data_agg, aes(x = value, col = colour, fill = colour)) +
        geom_histogram(alpha = 0.3, position = 'identity')+
        #geom_text(data = p_value_agg, aes(x=x, y=y, label = p_value)) +
        scale_fill_colorblind2() +
        scale_color_colorblind2() +
        theme(legend.position = 'none')

      if (isTRUE(log)) {
        p2 <- p2 + scale_x_log10()
      }


      if (isTRUE(ks_test)) p2 <- p2 +   geom_text(data = p_value_agg,
                                                  aes(x=x, y=y, label = p_value))

      if (return == 'both') return(p2 / p1)
      else if (return == 'agg') return(p2)
      else if (return == 'shares') return(p1)
      else stop('not implemented yet. set `return` to either "both", "agg", or "shares" ')
      #cowplot::plot_grid(p2, p1, ncol = 1)

    }



  } else if (type == 'xy') {
    sample_orig[, rank := frank(value), by = Var2]
    setorder(sample_orig, Var2, rank)
    sample_new[, rank := frank(value), by = Var2]
    setorder(sample_new, Var2, rank)

    data_ranked <- merge(sample_orig, sample_new, by = c('Var2', 'rank'),
                         suffixes = c("_orig", '_new'))

    # plot
    ggplot(data_ranked, aes(x = value_orig, y = value_new)) +
      geom_point(col = 'grey30') +
      facet_wrap(~Var2, scales = 'free') +
      geom_abline(slope = 1, intercept = 0)

  } else {
    stop('type not implemented yet')
  }
}

#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
compare_samples2 <- function(x, y) {
  summary_x <- apply(x,2, summary2)
  summary_y <- apply(y,2, summary2)
  summary_total_x <- summary2(rowSums(x))
  summary_total_y <- summary2(rowSums(y))

  reldif <- (summary_y - summary_x) / summary_x
  corx <- cor(x)
  cory <- cor(y)
  message("Mean and SD (disaggs):")
  print(sideBySide(list('Original sample:' = summary_x,
                        "New sample:"= summary_y,
                        "Reldif:" = reldif)))

  message("Mean and SD (agg):")
  print(sideBySide(list('Original sample:' = summary_total_x,
                        "New sample:"= summary_total_y,
                        "Reldif:" = (summary_total_y - summary_total_x) / summary_total_x)))

  message("Correlations")
  print(sideBySide(list(corx, cory)))
  cat('Relative difference:\n')
  print((cory - corx) / corx)
}



#' Title
#'
#' @param sample a matrix of samples. ncol(sample) = number of variables, nrow(sample) = number of samples
#' @param dist the dist to sample data from, currently only 'lnorm'
#' @param method the method used by `fitdistrplus::fitdist` (only relvant if skew of at least one variable is > max_skew)
#' @param max_skew the maximum skew above which distribution parameter are fitted using `fitdistrplus::fitdist`
#' @param type only if dist == 'truncnorm': which sampler to use (package name),
#' default: 'tmvtnorm'. for more efficient (but slighly less accurate results) choose 'tmvnsim'
#' @param algorithm only if dist == 'truncnorm': which algorithm to use (default: gibbs, see ?rtmvnorm)
#' @param fit only if dist == 'lnorm': should a lornorm dist be fitted to skewed samples using `fitdist`
#'
#' @return
#' @export
#'
#' @examples
reconstruct_disagg_sample <- function(sample, dist = 'lnorm',
                                      method = 'mme',
                                      type = NULL,
                                      algorithm = 'gibbs',
                                      fit = FALSE,
                                      max_skew = 2) {



  # calculate summary statistics for sample
  N <- nrow(sample)

  if (dist == 'lnorm' && (type == 'rlnorm.rplus' | type == 'univariate'))  {
    sample <- log(sample)
  }

  mean_orig <- apply(sample, 2, mean)
  sd_orig <- apply(sample, 2, sd)
  var_orig <- apply(sample, 2, var)
  #covmat_orig <- cov.shrink(sample)
  covmat_orig <- cov(sample)
  cormat_orig <- cor(sample)


  if (isTRUE(fit)) {
    skewness <- sapply(1:ncol(sample), function(x) skewness(sample[,x]))
    if (dist == 'lnorm' & type != 'mvlognormal' & length(skewness[skewness > max_skew]) > 0) {
      # at least one disagg variable has a high enough skew
      model_estimates <- sapply((1:ncol(sample))[skewness > max_skew], function(x) {
        result <- try(fitdistrplus::fitdist(sample[, x], distr = dist, method = method),
                      silent = TRUE)
        # TODO: which method delivers best fit??

        if (inherits(result, "try-error")) {
          warning(paste("Error in fitdist for column", x, "- using default parameter instead"))
          # Replace 'alternative_function' with your fallback function
          alternative_result <- c(
            mean_orig[x],
            diag(covmat_orig)[x]
          )
          return(alternative_result)
        } else {
          return(result$estimate)
        }
      })


      # replace diag of cov mat
      diag(covmat_orig)[skewness > max_skew] <- model_estimates[2,]
      mean_orig[skewness > max_skew] <- model_estimates[1,]

      # covmat_orig[!is.finite(covmat_orig)]
      # try(rlnorm.rplus(N, meanlog = mean_orig[skewness > max_skew],
      #              covmat_orig[skewness > max_skew, skewness > max_skew]))
      # try(chol.default(covmat_orig, pivot = TRUE))

    }
  }

  # sample
  if (dist == 'lnorm') {
    if (is.null(type) || type == 'rlnorm.rplus') {
      sample_new <- rlnorm.rplus(N, meanlog = mean_orig, covmat_orig)
      # TODO: performance comparison with other multi variate log normal sampler
      class(sample_new) <- 'matrix'
    } else if (type == 'mvlognormal') {
      sample_new <- mvlognormal(N, Mu = mean_orig, Sigma = (var_orig),
                                R = cor(log(sample)))
    } else if (type == 'univariate') {
      sample_new <- mapply(rlnorm, meanlog = mean_orig, sdlog = sd_orig, n = N)
      # sample_new <- rlnorm(N, meanlog = mean_orig, sdlog = sd_orig)
    } else if (type == 'univariate2') {
      sample_new <- mapply(rlnorm2, mean = mean_orig, sd = sd_orig, n = N)
      # sample_new <- rlnorm(N, meanlog = mean_orig, sdlog = sd_orig)
    }
  } else if (dist == 'truncnorm') {
    if (is.null(type) || type == 'tmvtnorm') {
      sample_new <- tmvtnorm::rtmvnorm(N, mean = mean_orig, sigma = covmat_orig,
                                       lower = rep(0, length = length(mean_orig)),
                                       upper = rep(Inf, length = length(mean_orig)),
                                       algorithm = algorithm)

    } else if (type == 'tmvnsim') {
      # fast
      sample_new <- tmvnsim::tmvnsim(N, k = length(mean_orig), means =  mean_orig,
                                     sigma = covmat_orig,
                                     lower = rep(0, length = length(mean_orig)),
                                     upper = rep(Inf, length = length(mean_orig)))$samp
    }
  } else if (dist == 'gamma'){
    sample_new <- rgamma2_vec(N, mean = mean_orig, sd = sd_orig)

  } else if (dist == 'mvgamma') {
    sample_new <- rmvgamma2(N, mean = mean_orig, sd = sd_orig, corr = cormat_orig)


  } else {
    stop('not implemtned yet')
  }
  return(sample_new)
}


rlnorm.rplus <- function (n, meanlog, varlog)
{
  D <- NCOL(oneOrDataset(meanlog))
  erg <- rplus(perturbe.aplus(exp(matrix(rnorm(n * length(meanlog)),
                                         ncol = D) %*% chol(varlog, pivot = TRUE)), exp(meanlog)))
  colnames(erg) <- colnames(oneOrDataset(mean))
  erg
}



plot1_full <- function(dt_agg, dt_shares,
                       colour_headers = 'grey70',
                       strip_text_size= 8, labels = TRUE) {

  if(isTRUE(labels)) {
    # create labels for plots
    dt_agg[dist2 == 'unif', dist := paste0('Unif(', min, ',', max, ')')]
    dt_agg[dist2 == 'exp', dist := paste0('Exp(1/',mean, ')')]
    dt_agg[dist2 == 'norm', dist := paste0('Norm(', mean, ',', sd, ')')]
    dt_agg[dist2 == 'truncnorm', dist := paste0('Truncnorm(', mean, ',', sd,',',min, ')')]


    dt_shares[, dist_disagg := dist_disagg2]
    dt_shares[dist_disagg2 != 'Dir_Generalised',
              dist_disagg := sapply(shares, function(x) {
                paste0('Dir(', paste(round(x,1), collapse = ','), ')')
              })]

    dt_shares[dist_disagg2 == 'Dir_Generalised',
              dist_disagg := mapply(function(x,y) {
                paste0('Dirg(', paste(x, collapse = ','),
                       ';\n',paste(y, collapse = ','),
                       ')')
              }, x = shares, y = sds)]

    dt_shares[dist_disagg2 == 'Dir_Nested',
              dist_disagg := mapply(function(x,y) {
                paste0('Dirn(', paste(x, collapse = ','),
                       ';\n',paste(y, collapse = ','),
                       ')')
              }, x = shares, y = sds)]

  } else {
    dt_agg[, dist := dist2]
    dt_shares[, dist_disagg := dist_disagg2]

  }

  # sample distributions
  sample_agg <- dt_sample_aggregates(dt_agg)
  sample_shares <- dt_sample_shares(dt_shares)

  # combine both
  sample_combined <- as.data.table(expand_grid(sample_agg, sample_shares))


  sample_combined[,
                  sample_disagg := pmap(list(sample_agg, sample_shares),
                                        function(x, y) x*y)]


  # Unnest samples
  agg <- unnest(sample_agg[, .(dist, sample_agg)],
                "sample_agg") %>%
    as.data.table
  shares <- unnest(sample_shares[, .(dist_disagg, sample_shares)],
                   "sample_shares") %>%
    as.data.table %>%
    melt(id.vars = c('dist_disagg'))
  disagg <- unnest(sample_combined[, .(dist, dist_disagg, sample_disagg)],
                   "sample_disagg") %>%
    as.data.table %>%
    melt(id.vars = c('dist', 'dist_disagg'))


  # Plot top row
  p_top <- ggplot(agg, aes(x = sample_agg)) +
    geom_histogram(alpha = 0.4, col = 'black', fill= 'black') +
    facet_wrap(~dist, nrow = 1, scales = 'free') +
    theme(axis.title = element_blank(),
          panel.background = element_rect(fill = 'grey80'),
          plot.background = element_rect(fill = colour_headers),
          strip.text = element_text(size = strip_text_size)
    )


  # Plot right column
  p_right <- ggplot(shares, aes(x = value, fill = variable, col = variable)) +
    geom_histogram(position = 'identity', alpha = 0.2) +
    facet_wrap(~dist_disagg, ncol = 1, scales = 'free', strip.position = 'right') +
    theme(axis.title.x = element_blank()) +
    scale_fill_colorblind7() +
    scale_color_colorblind7() +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          panel.background = element_rect(fill = 'grey80'),
          plot.background = element_rect(fill = colour_headers),
          strip.text = element_text(size = strip_text_size)

    )

  # plot inner part
  p_center <- ggplot(disagg, aes(x = value, fill = variable, col = variable)) +
    geom_histogram(position = 'identity', alpha = 0.2) +
    # stat_cor(aes(y = value, col = variable), label.x = 0, label.y = 600,
    #          cor.coef.name = 'rho') +
    # geom_text(data = result, aes(label = correlation_matrix),
    #           x = 20, y = 600, size = 3,
    #           parse = FALSE, inherit.aes = FALSE) +
    facet_grid(dist_disagg~dist, scales = 'free') +
    xlab('variable value') +
    scale_fill_colorblind7() +
    scale_color_colorblind7() +
    theme(legend.position = 'none',
          #panel.background = element_rect(fill = 'grey80'),
          #axis.title.x = element_blank(),
          strip.text = element_blank())


  # plot top right
  p_topright <- ggplot(data.frame(x = c(0, 1), y = c(0, 1)), aes(x, y)) +
    geom_blank() +
    theme(panel.background = element_rect(fill = colour_headers),
          panel.border = element_rect(fill = NA, colour= colour_headers),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_blank(),
          plot.background = element_rect(fill = colour_headers, colour = 'white')) +
    geom_abline(intercept = 0, slope = 1, color = "white", linewidth = 2) +
    annotate("text", x = 0.25, y = 0.75, label = "aggregate", size = 5,
             color = "black", hjust = 0) +
    annotate("text", x = 0.75, y = 0.25, label = "shares", size = 5,
             color = "black") +
    geom_segment(aes(x = 0.2, y = 0.75, xend = 0.0, yend = 0.75),
                 arrow = arrow(angle = 10, type = "closed", length = unit(0.2, "inches")),
                 color = "black") +
    geom_segment(aes(x = 0.75, y = 0.2, xend = 0.75, yend = 0.0),
                 arrow = arrow(angle = 10, type = "closed", length = unit(0.2, "inches")),
                 color = "black")

  # full plot
  p_top + wrap_elements(full = p_topright, clip = FALSE) + #(plot_spacer() + theme(plot.background = element_rect(fill = 'grey70'))) +
    p_center + p_right +
    plot_layout(nrow = 2, ncol = 2, widths = c(4,1), heights = c(1,4))

}


plot2_full <- function(dt_agg, dt_shares) {


  dt_agg[, dist2 := factor(dist2, levels = c('norm', 'exp', 'unif'))]
  # dt_shares[, dist_disagg2 := factor(dist_disagg2,
  #                                      levels = c('dir1', 'dir_maxent', 'gdir'))]


  dt_agg[dist2 == 'unif', dist := paste0('Unif(', min, ',', max, ')')]
  dt_agg[dist2 == 'exp', dist := paste0('Exp(1/',mean, ')')]
  dt_agg[dist2 == 'norm', dist := paste0('Norm(', mean, ',SD)')]
  #dt_agg[dist2 == 'truncnorm', dist := paste0('Truncnorm(', mean, ',', sd,',',a, ')')]

  #dt_shares[, dist_disagg := dist_disagg2 %>% as.factor]
  #dt_shares[, dist_disagg := dist_disagg2 %>% as.factor]

  dt_shares[is.na(sds),
            dist_disagg := sapply(shares, function(x) {
              paste0('Dir(', paste(round(x, 2), collapse = ','), ')')
            })]

  dt_shares[!is.na(sds),
            dist_disagg := mapply(function(x,y) {
              paste0('Dirg(', paste(x, collapse = ','),
                     ';CV)')
            }, x = shares, y = sds)]


  # dt_shares[!(dist_disagg2 %in% c('gdir','gdir1', 'gdir2', 'gdir3')),
  #             dist_disagg := sapply(alpha, function(x) {
  #               paste0('Dir(', paste(x, collapse = ','), ')')
  #             })]
  #
  # dt_shares[dist_disagg2 %in% c('gdir','gdir1', 'gdir2', 'gdir3'),
  #             dist_disagg := mapply(function(x,y) {
  #               paste0(dist_disagg2, '(', paste(x, collapse = ','),
  #                      ';beta)')
  #             }, x = alpha, y = beta) %>% factor]




  # sample distributions
  sample_agg2 <- dt_sample_aggregates(dt_agg)
  sample_shares2 <- dt_sample_shares(dt_shares)


  # combine both
  sample_combined2 <- as.data.table(expand_grid(sample_agg2, sample_shares2))
  sample_combined2[,
                   sample_disagg := pmap(list(sample_agg, sample_shares),
                                         function(x, y) x*y)]
  sample_combined2[, id := 1:.N]

  # Calculate correlation coefficienct rho
  pars <- sample_combined2
  #pars <- pars_sample
  pars[, rho := lapply(sample_disagg, function(x) {
    rho <- cor(x)
    return(rho[lower.tri(rho)])
  })]
  pars <- pars[dist != 'truncnorm']

  # unnest data.table
  data <- as.data.table(unnest(pars[, .(id, mean, sd, min, max, shares,
                                        sds, dist,
                                        dist2, dist_disagg2,
                                        dist_disagg, rho)],
                               'rho'))

  # create labels where necessary
  data[, xlabel := '']
  data[, ylabel := '']
  # data[dist2 == 'unif', `:=`(xvar = b - a,
  #                           xlabel = 'b-a')]

  data[!is.na(sd), `:=`(xvar = sd,
                        xlabel = 'SD')]
  # data[dist2 == 'exp', `:=`(xvar = 1/mean,
  #                          xlabel = '1/mean')]

  data[!is.na(sds), `:=`(yvar = sapply(sds, function(x) x[1]),
                         ylabel = 'CV')]




  #plot_correlation_bubble <- function(data, x = 'x', y = 'y')

  data[, abs_value := abs(rho)]

  #data[, x := as.factor(sd)]
  #data[, y := as.factor(beta2)]

  # set the order
  setorder(data, id, xvar, yvar, -abs_value)
  data[, id2 := rep(1:3, nrow(data) / 3)]


  data[, x := as.numeric(as.factor(xvar)), by = .(dist, dist_disagg)]
  data[, y := as.numeric(as.factor(yvar)), by = .(dist, dist_disagg)]

  unique_x <- unique(na.omit(data$x))
  unique_y <- unique(na.omit(data$y))

  data[is.na(x), x := mean(unique_x)]
  data[is.na(y), y := mean(unique_y)]

  # Scaling the sizes of the bubbles
  max_size <- 0.5
  min_size <- 0.01
  data$size <- calculate_size(data$rho)
  data$radius <- data$size / 2
  data[id2 == 1, xpos := as.double(x)]
  data[id2 == 1, ypos := y + radius]
  data[id2 == 2, xpos := x - (radius / sqrt(2))]
  data[id2 == 2, ypos := y - (radius / sqrt(2))]
  data[id2 == 3, xpos := x + (radius / sqrt(2))]
  data[id2 == 3, ypos := y - (radius / sqrt(2))]

  # assigning pane number (1 to 9), row-wise
  panel_number <- as.data.table(expand.grid(dist = unique(data$dist),
                                            dist_disagg = unique(data$dist_disagg)))
  panel_number[, pane_id := 1:.N]
  data <- merge(data, panel_number, by = c('dist', 'dist_disagg'),
                all = TRUE, sort = FALSE)


  # split by pane number
  setorder(data, pane_id)
  data2 <- split(data, by = c('pane_id'), flatten = TRUE)
  #panel_names <- str_split(names(data2), '\\.')


  # data2 <- split(data, by = c('dist', 'dist_disagg'), flatten = TRUE)
  # panel_names <- str_split(names(data2), '\\.')

  plist <- vector('list', length(data2))
  #plist <- setNames(plist, names(data2))
  offset <- 0.5


  theme_set(theme_bw())
  for (i in 1:length(data2)) {
    #cat(i, '')
    #name <- str_split(names(plist)[i], '\\.') %>% unlist

    plist[[i]] <- ggplot(data2[[i]], aes(x0 = xpos, y0 = ypos,
                                         r = radius,fill = rho)) +
      geom_circle(alpha = 1, color = NA) +
      scico::scale_fill_scico(palette = "vik", midpoint = 0, direction = -1L,
                              limit = c(-1, 1)) +
      #theme_bw() +
      xlab(unique(data2[[i]]$xlabel)) +
      ylab(unique(data2[[i]]$ylabel)) +
      scale_x_continuous(breaks = unique(data2[[i]]$x), labels = data2[[i]]$xvar %>% unique,
                         limits = c(min(unique_x) - offset, max(unique_x) + offset)) +
      scale_y_continuous(breaks = unique(data2[[i]]$y), labels = data2[[i]]$yvar %>% unique,
                         limits = c(min(unique_y) - offset, max(unique_y) + offset)) +
      coord_fixed(ratio = 1) +
      theme(legend.position = 'none',
            #axis.text = element_blank(),
            #axis.ticks = element_blank(),
            #axis.title = element_blank(),
            #aspect.ratio = 1,
            plot.margin = margin(1,1,1,1),
            panel.border = element_rect(color = 'grey30'),
            strip.background  = element_rect(fill = 'grey90', colour = NA),
            panel.grid.major =  element_blank()
      )
    #plist[[i]]
    if (i %in% 1:2) {
      # top row
      plist[[i]] <- plist[[i]] +
        facet_wrap(~dist)
    }
    if (i %in% c(6,9)) {
      # right column
      plist[[i]] <- plist[[i]] +
        facet_wrap(~dist_disagg, strip.position = 'right')
    }
    if (i == 3) {
      # topright
      plist[[i]] <- plist[[i]] +
        facet_grid(dist_disagg ~ dist)
    }

    if (!(i %in% c(7:9))) {
      # bottom
      plist[[i]] <- plist[[i]] +
        theme(axis.text.x = element_blank(),
              #axis.ticks.x = element_blank(),
              axis.title.x = element_blank() )
    }
    if (!(i %in% c(1,4,7))) {
      # left
      plist[[i]] <- plist[[i]] +
        theme(axis.text.y = element_blank(),
              #axis.ticks.y = element_blank(),
              axis.title.y = element_blank() )
    }
    if (unique(data2[[i]]$ylabel) == '') {
      plist[[i]] <- plist[[i]] +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    }
    if (unique(data2[[i]]$xlabel) == '') {
      plist[[i]] <- plist[[i]] +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

    }
    # plist[[i]] <- plist[[i]] + theme_dark()
  }
  # legend

  legend <- get_legend(ggplot(data = data.table(rho = c(-1, 1)),
                              aes(fill = rho)) +
                         scico::scale_fill_scico(palette = "vik", midpoint = 0, direction = -1L) +
                         geom_point(x = 1, y = 1) +
                         theme(legend.direction = 'vertical',
                               plot.margin = margin(0,0,0,0))
  )
  # legend <- as_ggplot(legend)
  #library(patchwork)
  p <- wrap_plots(plist, ncol = 3, nrow = 3, byrow = TRUE)
  wrap_plots(p, wrap_elements(full = legend), ncol = 2, widths = c(10,1))

  # p + inset_element(legend, left = 1, bottom = 1, right = 1.5, top = 1,
  #                   align_to = 'full', on_top = FALSE)
  #
}

plot3full <- function(dt_agg, dt_shares,
                      sample_dist = 'lnorm',
                      sample_type = 'mvlognormal',
                      type = 'hist',
                      return = 'both') {

  # create labels
  dt_agg[dist2 == 'unif', dist := paste0('Unif(', a, ',', b, ')')]
  dt_agg[dist2 == 'exp', dist := paste0('Exp(1/',mean, ')')]
  dt_agg[dist2 == 'norm', dist := paste0('Norm(', mean, ',', sd, ')')]
  dt_agg[dist2 == 'truncnorm', dist := paste0('Truncnorm(', mean, ',', sd,',',a, ')')]

  dt_shares[dist_disagg2 != 'Dir_Generalised'
            , dist_disagg := sapply(alpha, function(x) {
              paste0('Dir(', paste(round(x,2), collapse = ','), ')')
            })]

  dt_shares[dist_disagg2 == 'Dir_Generalised',
            dist_disagg := mapply(function(x,y) {
              paste0('Dirg(', paste(x, collapse = ','),
                     ';\n',paste(y, collapse = ','),
                     ')')
            }, x = alpha, y = beta)]

  # sample both, agg and shares and combine both to get sample of the disagg
  sample_agg <- dt_sample_aggregates(dt_agg)
  sample_shares <- dt_sample_shares(dt_shares)

  sample_combined <- as.data.table(expand_grid(sample_agg, sample_shares))
  sample_combined[, sample_disagg := pmap(list(sample_agg, sample_shares),
                                          function(x, y) x*y)]


  # plot
  sample_combined[, id := 1:.N]

  # resample
  sample_combined[, sample_lnorm := (lapply(sample_disagg,reconstruct_disagg_sample,
                                            dist = sample_dist,
                                            type = sample_type))] #mvlognormal
  if (type == 'hist') {
    # histogram
    sample_combined[, plot := Map(f = compare_samples,
                                  sample_orig = sample_disagg,
                                  sample_new = sample_lnorm,
                                  ks_test = FALSE,
                                  print = F,
                                  return = return,
                                  #log = TRUE,
                                  type = 'hist')]

    sample_combined[, plot := lapply(plot, function(x) {
      wrap_elements(
        # x & theme(plot.background = element_rect(fill = 'grey90',
        #                                          colour = 'grey90'))
        full = x + plot_annotation(theme = theme_border) &
          xlab(NULL) & ylab(NULL),
        clip = FALSE
      )
      #x & theme(plot.background = element_rect(colour = 'red'))

      #x + plot_annotation(theme = theme_border)
    }) ]

  } else if (type == 'ecdf') {
    # ecdf
    sample_combined[, plot := Map(f = compare_samples,
                                  sample_orig = sample_disagg,
                                  sample_new = sample_lnorm,
                                  ks_test = FALSE,
                                  print = FALSE,
                                  return = return,
                                  type = 'ks')]

    sample_combined[ ,#& even(id) == TRUE,
                     plot := lapply(plot, function(x) {
                       wrap_elements(
                         # x & theme(plot.background = element_rect(fill = 'grey90',
                         #                                          colour = 'grey90'))
                         full = x + plot_annotation(theme = theme_border) &
                           xlab(NULL) & ylab(NULL) & theme(legend.position = 'none'),
                         clip = FALSE
                       )
                       #x & theme(plot.background = element_rect(colour = 'red'))

                       #x + plot_annotation(theme = theme_border)
                     }) ]

  } else {
    stop('type not implemented')
  }

  p_center <- wrap_plots(sample_combined$plot,
                         ncol = 3, nrow = 3, guides = 'collect', byrow = FALSE)
  p_top <- gg_strip2(unique(sample_combined$dist))
  p_right <- gg_strip2(unique(sample_combined$dist_disagg),
                       direction = 'y')


  # wrap plot together
  wrap_plots(plot_spacer(), p_top, plot_spacer(),
             grid::textGrob('count', rot = 90), p_center, (p_right),
             plot_spacer(), grid::textGrob('value'), plot_spacer(),
             ncol = 3,
             widths = c(0.5, 15,1),
             heights = c(1,15,0.5)
  ) # TODO: add legend

}




# Fill
scale_fill_colorblind7 = function(.ColorList = 2L:8L, ...){
  scale_fill_discrete(..., type = colorblind_pal()(8)[.ColorList])
}

# Color
scale_color_colorblind7 = function(.ColorList = 2L:8L, ...){
  scale_color_discrete(..., type = colorblind_pal()(8)[.ColorList])
}

# Fill
scale_fill_colorblind2 = function(.ColorList = c(1L, 2L), ...){
  scale_fill_discrete(..., type = colorblind_pal()(8)[.ColorList])
}

# Color
scale_color_colorblind2 = function(.ColorList = c(1L, 2L), ...){
  scale_color_discrete(..., type = colorblind_pal()(8)[.ColorList])
}

# Fill
scale_fill_colorblind3 = function(.ColorList = c(1L, 3L:8L), ...){
  scale_fill_discrete(..., type = colorblind_pal()(8)[.ColorList])
}

# Color
scale_color_colorblind3 = function(.ColorList = c(1L, 3L:8L), ...){
  scale_color_discrete(..., type = colorblind_pal()(8)[.ColorList])
}


plot_full <- function(sample_agg, sample_shares, sample_disagg) {
  sample_agg <- as.data.table(sample_agg)
  sample_shares <- as.data.table(sample_shares)
  sample_disagg <- as.data.table(sample_disagg)

  setnames(sample_agg, 'agg')
  # save default plot settings
  par_default <- par()

  # plot histograms
  p_agg <- ggplot(sample_agg, aes(x = agg)) +
    geom_histogram(alpha = 0.4)

  p_shares <- ggplot(melt(sample_shares) %>% suppressWarnings(),
                     aes(x = value, col = variable, fill = variable)) +
    geom_histogram(position = 'identity', alpha = 0.4) +
    scale_fill_colorblind7() +
    scale_color_colorblind7() +
    theme(legend.position = 'none')

  p_disagg <- ggplot(melt(sample_disagg) %>% suppressWarnings(),
                     aes(x = value, col = variable, fill = variable)) +
    geom_histogram(position = 'identity', alpha = 0.4) +
    scale_fill_colorblind7() +
    scale_color_colorblind7() +
    theme(legend.position = 'none')

  dev.control('enable') # enable display list
  # plot ternary plot for shares
  par(mar = rep(0, 4), oma = c(0, 0, 0, 0), bg = rgb(0,0,0,0))
  {
    TernaryPlot(lab.col = my_cols[2:4],
                ticks.col = my_cols[2:4],
                grid.col = my_cols[2:4],
                axis.labels = seq(0, 1, by = .1))
    TernaryPoints(sample_shares, pch = 16, cex = 0.5,
                  col= grey(level = 0.3, alpha = 0.4))
  }
  p_tern <- recordPlot()

  # plot correaltion matrix
  par(oma = c(0,0,0, 2), bg = rgb(0,0,0,0))
  cor <- cor(sample_disagg)
  corrplot(cor, diag = FALSE, tl.pos = 'd', tl.col = my_cols[2:4],
           tl.cex = 2)
  p_cor <- recordPlot()

  dev.off()

  # combine plots
  p_shares_comb <- plot_grid(p_shares, p_tern, ncol = 1, align = 'r')
  p_disagg_comb <- plot_grid(p_disagg, p_cor, ncol = 1)

  #par( mar = rep(3, 4), oma = rep(3, 4))
  par(mar = c(5, 4, 4, 2) + 0.1, oma = c(3,3,3,3), bg = 'white')
  return(plot_grid(p_agg, p_shares_comb, p_disagg_comb,
                   nrow = 1, labels = c('Agg', 'Shares', 'Disagg'),
                   vjust = 0))

}



ternary_hists <- function(sample, oma = 5, mar = 4,
                          bins = seq(0,1,0.025),
                          hist_height = 1) {
  # calculate hist stats
  hist_stats <- lapply(1:3, function(x) {
    p <- hist(sample[,x], plot = FALSE)
    data.frame('counts' = max(p$counts),
               'nbins' = length(p$mids))
  } )
  hist_stats <- do.call(rbind, hist_stats)

  max_count <- max(hist_stats[,'counts'])
  max_nbin <- max(hist_stats[,'nbins'])

  if (length(bins) == 1 && bins == 'max_nbin') bins = max_nbin

  # plot histograms and save to tempdir
  tempdir <- tempdir(check = TRUE)
  p_hists <- lapply(1:3, function(x) {
    png(file.path(tempdir, paste0('temp', x, '.png')), bg = 'grey')
    par(oma = c(0, 0, 0, 0), bg = rgb(0,0,0,0),
        mar = c(mar, 0,0,0), xpd = NA)
    hist(sample[,x], axes = FALSE, main = NULL,
         xlab = NULL, ylab = NULL,
         xlim = c(0,1),
         breaks = bins,
         ylim = c(0, max_count *(1/hist_height)),
         border = my_cols[x+1],
         col = adjustcolor(my_cols[x+1], alpha.f = 0.5))
    abline(h = 0, col = 'white', lwd = 1.3)
    # axis(1, seq(0, 1, 0.1), labels=FALSE, col = 'white',
    #      col.ticks = my_cols[x+1],
    #      tck=-0.5)



    #abline(v = seq(0, 1, .1), lty = 2, col = my_cols[x+1])
    dev.off()

  } )

  # read in hists againas image
  hists <- lapply(list.files(tempdir, pattern = '^temp.*\\.png$', full.names = TRUE),
                  image_read)

  # rotate hists
  hists[[1]] <- image_rotate(hists[[1]], -60)
  hists[[2]] <- image_rotate(hists[[2]], 60)
  hists[[3]] <- image_rotate(hists[[3]], 180)

  # make background transparent
  hists <- lapply(hists, image_transparent, color = 'white')

  # plot ternarny
  par(mar = rep(oma, 4))
  TernaryPlot(lab.col = my_cols[2:4],
              ticks.col = my_cols[2:4],
              grid.col = my_cols[2:4],
              axis.labels = seq(0, 1, by = .1))
  TernaryPoints(sample, pch = 16, cex = 0.5,
                col= grey(level = 0.3, alpha = 0.4))

  # add histograms
  point <- 0.8660254 + 1/sqrt(3)
  rasterImage(hists[[1]], -point,0,0,point, xpd=TRUE)
  rasterImage(hists[[2]], 0,0,point,point, xpd=TRUE)
  rasterImage(hists[[3]], -0.5,-(480/658)*point,0.5,0, xpd=TRUE)
}





is_truncnorm <- function(mean, sd, a, b) {
  return(is.numeric(mean) & is.numeric(sd) & (is.finite(a) | is.finite(b)))
}

is_norm <- function(mean, sd, a, b) {
  return(is.numeric(mean) & is.numeric(sd) & is.infinite(a) & is.infinite(b))
}
is_exp <- function(mean, sd, a, b) {
  return(is.numeric(mean) & !is.numeric(sd) & a == 0)
}
is_unif <- function(mean,sd, a, b) {
  return(!is.numeric(mean) & is.finite(a) & is.finite(b))
}

choose_dist_agg <- function(mean, sd, a, b) {
  if (is.finite(a) & is.finite(b) & a >= b) {
    warning('a must be < b')
    return(NA)
  }
  if (is.finite(sd) & sd <=0) {
    warning('sd must be positive')
    return(NA)
  }
  if (is.finite(mean) & is.finite(a) & mean <= a) {
    warning('mean must be > a')
    return(NA)
  }
  if (is.finite(mean) & is.finite(b) & mean >= b) {
    warning('mean must be < b')
    return(NA)
  }

  if (is.finite(mean)) {
    if (is.finite(sd)) {
      if (is.finite(a) | is.finite(b)) {
        return('truncnorm')
      } else {
        return('norm')
      }
    } else if (isTRUE(all.equal(a, 0)) & !is.finite(b)) {
      return('exp')
    } else {
      return(NA)
    }
  } else if (is.finite(a) & is.finite(b) & !is.finite(sd)) {
    return('unif')
  } else {
    NA
  }
}


sample_par_comb <- function(x, plot = FALSE) {
  x[, dist := pmap_chr(list(mean, sd, a, b), choose_dist_agg)]
  x <- x[!is.na(dist)]


  x[, dist_disagg := pmap_chr(list(alpha, beta), function(alpha, beta) {
    if(length(beta) == 1 && is.na(beta)) beta <- NULL
    choose_dist_disagg(alpha, beta)
  } )]
  x <- x[!(sapply(alpha, var) == 0 & !is.na(beta))] #exclude combination uninformative alpha and beta

  # sample aggregate
  x[, sample_agg := pmap(list(mean, sd, a, b), function(mean, sd, a, b) {
    if(is.na(sd)) sd <- NULL
    if(is.na(mean)) mean <- NULL
    tryCatch(
      ragg(N, mean, sd, a, b),
      error = function(e) e
    )
  })]

  # sample shares
  x[, sample_shares := pmap(list(alpha, beta), function(alpha, beta) {
    if(length(beta) == 1 && is.na(beta)) beta <- NULL
    tryCatch(
      rshares(N, alpha, beta),
      error = function(e) e
    )
  })]

  # combine both
  x[!sapply(sample_agg, function(x) inherits(x, 'error')),
    sample_disagg := pmap(list(sample_agg, sample_shares), function(x, y) x*y)]

  # check which combinations still have to be implemented
  errors <- x[sapply(sample_agg, function(x) inherits(x, 'error'))]
  if (nrow(errors) > 0) {
    warning('those combinations could be sampled:')
    print(errors)
  }

  # create ids
  x[, id := 1:.N]
  setcolorder(x, 'id')

  # plot
  if (isTRUE(plot)) {
    x[, plot := pmap(list(sample_agg, sample_shares, sample_disagg),
                     plot_full)]

  }
  return(x[])
}


dt_sample_aggregates <- function(x, N = 10000) {
  #x[, dist := pmap_chr(list(mean, sd, a, b), choose_dist_agg)]
  #x <- x[!is.na(dist)]

  # sample aggregate
  x[, sample_agg := pmap(list(mean, sd, min, max), function(mean, sd, min, max) {
    if(is.na(sd)) sd <- NULL
    if(is.na(mean)) mean <- NULL
    tryCatch(
      ragg(N, mean, sd, min, max),
      error = function(e) e
    )
  })]
  return(x[])
}

dt_sample_shares <- function(x, N = 10000) {
  # x[, dist_disagg := pmap_chr(list(alpha, beta), function(alpha, beta) {
  #   if(length(beta) == 1 && is.na(beta)) beta <- NULL
  #   choose_dist_disagg(alpha, beta)
  # } )]
  # x <- x[!(sapply(alpha, var) == 0 & !is.na(beta))] #exclude combination uninformative alpha and beta


  # sample shares
  x[, sample_shares := pmap(list(shares, sds), function(shares, sds) {
    if(length(sds) == 1 && is.na(sds)) sds <- NULL
    tryCatch(
      rshares(N, shares, sds),
      error = function(e) e
    )
  })]
  return(x[])
}




# Function to calculate correlation matrix
calculate_correlation <- function(subset_data) {
  cor <- cor(subset_data[, .(sample_disagg.V1, sample_disagg.V2, sample_disagg.V3)])
  cor <- cor[lower.tri(cor)]
  cor <- round(cor, 2)
  text <- paste0('rho[1,2] = ', cor[1], '\n',
                 'rho[1,3] = ', cor[2], '\n',
                 'rho[2,3] = ', cor[3], '\n')

  return(text)
}

#' Title
#'
#' @param labels a vector of labels
#' @param label_size
#' @param label_col
#' @param background_fill
#' @param direction
#'
#' @return
#' @export
#'
#' @examples
gg_strip <- function(labels, label_size = 8, label_col = 'black',
                     background_fill = 'grey', direction = 'x') {
  if (direction == 'x') {
    nrow = 1
    ncol = NULL
    angle = 0
  } else if (direction == 'y') {
    nrow = NULL
    ncol = 1
    angle = 270
  }
  plots <- lapply(labels,
                  function(x) {
                    ggplot() +
                      # Adding a grey background
                      theme(panel.background = element_rect(fill = background_fill, colour = NA),
                            # Removing axis titles, text, and ticks
                            axis.title = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            # Removing panel grid, legend, and title
                            panel.grid = element_blank(),
                            legend.position = "none",
                            plot.title = element_blank()) +
                      # Adding the word 'foo' in the middle of the plot
                      annotate("text", x = 0.5, y = 0.5, label = x,
                               size = 8, color = label_col, angle = angle)

                  })

  wrap_plots(plots, nrow = nrow, ncol = ncol)

}

gg_strip2 <- function(labels, label_size = 12, label_col = 'black',
                      background_fill = 'grey', direction = 'x') {
  if (direction == 'x') {
    nrow = 1
    ncol = NULL
    angle = 0
  } else if (direction == 'y') {
    nrow = NULL
    ncol = 1
    angle = 270
  }

  grobs <- lapply(labels, function(x) {
    wrap_elements(full = grobTree( rectGrob(gp=gpar(fill=background_fill, col = 'white')),
                                   textGrob(x,
                                            gp=gpar(fontsize=label_size, col=label_col),
                                            rot = angle)),
                  clip = FALSE)

  })

  wrap_plots(grobs, nrow = nrow, ncol = ncol)
}


calculate_size <- function(x, max_size = 0.5, min_size = 0.01) {
  (abs(x) - min(abs(x))) / (max(abs(x)) - min(abs(x))) * (max_size - min_size) + min_size
}


#' Find Subclassifcations in a Hierarchical Structure
#'
#' This function searches for all direct subcategories of a given node in a hierarchical structure.
#' The hierarchy is defined by strings separated by pipe characters.
#' It returns only the immediate subcategories of the specified node, excluding the node itself and any deeper subcategories.
#'
#' @param all_codes A character vector containing the hierarchical codes.
#' @param node A character string representing the node for which direct subcategories are to be found.
#'
#' @return A character vector containing the direct subcategories of the given node.
#' @export
#'
#' @examples
#' all_codes <- c("total_for_category|liquid_fuels|gasoline", "total_for_category|liquid_fuels",
#'                "total_for_category|solid_fuels", "total_for_category", "total_for_category|biomass",
#'                "total_for_category|gaseous_fuels", "total_for_category|liquid_fuels|diesel_oil",
#'                "total_for_category|liquid_fuels|gas-diesel_oil")
#' node <- "total_for_category|liquid_fuels"
#' find_subclassifications(all_codes, node)
find_subclassifications <- function(all_codes, node) {
  # Escape any special characters in the node, including the pipe character
  escaped_node <- gsub("\\|", "\\\\|", node)

  # Create a pattern that matches only codes that are direct subcategories of the given node
  pattern <- paste0("^", escaped_node, "\\|[^|]+$")

  # Filter codes that match this pattern
  sub_codes <- all_codes[grepl(pattern, all_codes)]

  return(sub_codes)
}

#' Find Hierarchical Codes Below a Given Node
#'
#' This function searches through a list of hierarchical codes and returns
#' all codes that are hierarchically below a specified node. The specified node
#' itself is not included in the results.
#'
#' @param all_codes A vector of strings, each representing a hierarchical code.
#' @param node A string representing the node for which sub-nodes are to be found.
#'
#' @return A vector containing all hierarchical codes that are below the specified node.
#'
#' @examples
#' all_codes <- c("I", "I.1", "I.1.A", "I.1.A.1", "I.1.A.1.a", "I.1.B", "I.2", "I.3")
#' node <- "I.1.A"
#' find_hierarchical_codes_below(all_codes, node)
#'
#' @export
find_subcategories <- function(all_codes, node) {
  # Create a pattern that matches only codes that are below the given node
  pattern <- paste0("^", node, "\\..+$")

  # Filter codes that match this pattern
  sub_codes <- all_codes[grepl(pattern, all_codes)]

  return(sub_codes)
}

#' Find Hierarchical Codes Above a Given Node
#'
#' This function searches through a list of hierarchical codes and returns
#' all codes that are ancestors of (i.e., above) a specified node.
#' The specified node itself is not included in the results.
#'
#' @param all_codes A vector of strings, each representing a hierarchical code.
#' @param node A string representing the node for which ancestor nodes are to be found.
#'
#' @return A vector containing all hierarchical codes that are ancestors of the specified node.
#'
#' @examples
#' all_codes <- c("I", "I.1", "I.1.A", "I.1.A.1", "I.1.A.1.a", "I.1.B", "I.2", "I.3")
#' node <- "I.1.A.1.a"
#' find_hierarchical_codes_above(all_codes, node)
#'
#' @export
find_parent_categories <- function(all_codes, node) {
  # Split the node into its hierarchical components
  node_components <- unlist(strsplit(node, "\\."))

  # Generate all ancestor codes
  ancestor_codes <- sapply(seq_along(node_components), function(i) {
    paste(node_components[1:i], collapse = ".")
  })

  # Filter out the node itself and return ancestors that are in the all_codes list
  ancestors <- ancestor_codes[ancestor_codes != node]
  valid_ancestors <- ancestors[ancestors %in% all_codes]

  return(valid_ancestors)
}

#' Convert Sample Matrix to Data Table
#'
#' This function converts a sample matrix into a long-format data.table. It melts the matrix
#' into a long format and then creates a list-column in the resulting data.table, where each
#' list element consists of values from the matrix, named according to their corresponding
#' 'y' values.
#'
#' @param sample A matrix or data.frame that needs to be converted to a long-format data.table.
#' @param remove_non_finite A logical value indicating whether non-finite values (NA, NaN, Inf)
#'        should be removed from the matrix. Defaults to TRUE.
#'
#' @return A data.table where each row corresponds to a unique 'x' and includes a list-column
#'         'sample'. Each element of the list is named according to the corresponding 'y' value.
#'
#' @examples
#' # Assuming matrix 'mat' with some sample data
#' # mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
#' # result <- convert_sample_to_dt(mat)
#' # print(result)
#'
#' @importFrom data.table melt
#' @importFrom data.table Map
#' @export
convert_sample_to_dt <- function(sample, remove_non_finite = TRUE) {
  dt_long <- melt_matrix(sample, remove_zero = FALSE,
                         remove_non_finite = remove_non_finite)


  result <- dt_long[, .(sample = list(value)), by = y]
  # result <- dt_long[, .(sample = Map(setNames, list(value),
  #                                    list(as.character(x)))), by = y]

  return(result)
}



#' Extracts the variables from a function generated by `disaggR::generate_sampling_fun`
#'
#' @param fun
#'
#' @return
#' @export
#'
#' @examples
get_vars <- function(fun) {
  vars <- ls(envir = environment(fun))
  res <- lapply(vars, get, envir = environment(fun))
  names(res) <- vars
  return(res)
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
normalize_UNFCCC_classfications <- function(x) {
  x <- tolower(x)
  x <- gsub(' ', '_', x)
  x <- gsub('\\/', '-', x)
  return(x)
}

#' Get use or supply shares from SUTs (result of script 0a_prepare_exiobase_SUTs.R)
#'
#' Extract use/supply shares for
#' - (a) specific product(s)
#' - from (a) specific countr(ies)
#' - in specific industries
#' - in specific countries
#'
#'
#' @param x a supply or use table (result of script 0a_prepare_exiobase_SUTs.R)
#' @param products
#' @param industries
#' @param product_countries
#' @param industry_countries
#'
#' @return
#' @export
#'
#' @examples
# sut = use
# products = temp$EXIOBASE_products %>% unlist %>% unique
# industries = temp$EXIOBASE_code %>% unlist %>% unique %>% as.factor()
# industry_countries = c('FIN', 'DEU')
# product_countries <- 'all'

get_SUT_shares <- function(sut,
                           products = 'all',
                           industries = 'all',
                           product_countries = 'all',
                           industry_countries = 'all') {
  #sut <- copy(sut)
  #print(sut$product_code)
  #print(is.data.table(sut))
  if (length(products) == 1 && products == 'all') {
    products <- sut$product_code
  }
  if (length(industries) == 1 && industries == 'all') industries <- sut$industry_code
  if (length(industry_countries) == 1 && industry_countries == 'all') industry_countries <- sut$country_industry
  SUT_shares <- sut[
    product_code %in% products
    & industry_code %in% industries
    & country_industry %in% industry_countries
    , list(value = sum(value))
    , by = .(industry_code)
  ]
  SUT_shares[, country_total := sum(value)]
  SUT_shares[, share := value / country_total]


  if (nrow(SUT_shares) == 0) {
    # product not used in any of the industries --> equally share them among all
    # TODO: more elaborate
    SUT_shares <- data.table(industry_code = industries,
                             share = 1 / length(industries))
  }

  return(SUT_shares[, .(industry_code, share)])

}

#' Calculate Summary Statistics
#'
#' This function calculates various summary statistics based on a numeric vector column in a data.table.
#'
#' @param dt A data.table object with a column named 'sample' containing numeric vectors of length N.
#' @return The input data.table with additional columns for mean, median, standard deviation, coefficient of variation,
#' relative mean, and confidence intervals.
#' @import data.table
#' @import purrr
#' @keywords data manipulation
#' @export
calculate_summary_statistics <- function(dt) {
  dt[, mean := map_dbl(sample, mean, na.rm = TRUE)]
  dt[, median := map_dbl(sample, median, na.rm = TRUE)]
  dt[, sd := map_dbl(sample, sd, na.rm = TRUE)]
  dt[, cv := sd / mean]
  #dt[, mean_rel := mean / sum(mean, na.rm = TRUE), by = gas]
  dt[, CI97.5 := map_dbl(sample, quantile, probs = 0.975, na.rm = TRUE)]
  dt[, CI2.5 := map_dbl(sample, quantile, probs = 0.025, na.rm = TRUE)]
  dt[, sample_norm := mapply(function(x,y) x / y,
                             x = sample,
                             y = mean,
                             SIMPLIFY = FALSE)]
  return(dt[])
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
sum_samples <- function(x) {
  list(Reduce(`+`, x))
}

read_EB3_S_meta <- function(path) {
  colnames <- fread(file.path(path), nrows = 2,
                    drop = c(1), header = FALSE) %>%
    t %>%
    as.data.table %>%
    setnames(new = c('region', 'sector'))

  rownames <- fread(file.path(path), skip = 26,
                    select = c(1)) %>%
    setnames(new = c('category'))

  return(list(colnames = colnames, rownames = rownames))

}


as_sparse_matrix_list <- function(dt, i = 'row', j = 'col', sample = 'sample',
                                  nrow, ncol) {
  if (missing(nrow)) nrow = max(dt[[i]])
  if (missing(ncol)) ncol = max(dt[[j]])

  n <- length(dt[[sample]][[1]])

  matlist <- lapply(1:n, function(x) {
    sparseMatrix(dt[[i]], dt[[j]], x = sapply(dt[[sample]], '[', x),
                 dims = c(nrow, ncol))
  })
  return(matlist)
}



convert_sparse_matrix_to_dt <- function(x) {
  # from: https://stackoverflow.com/questions/21099612/extract-i-and-j-from-a-sparse-matrix
  data.table(
    i = as.integer(x@i + 1),
    # Cols (a little tricky..)
    j = as.integer(findInterval(seq(x@x)-1,x@p[-1]) + 1),
    # Values
    value = x@x
  )
}

convert_dt_to_sparse_matrix <- function(dt, i = 'i', j = 'j', x = 'value',
                                        nrow, ncol) {
  if (missing(nrow)) nrow = max(dt[[i]])
  if (missing(ncol)) ncol = max(dt[[j]])

  sparseMatrix(dt[[i]], dt[[j]], x = dt[[x]],
               dims = c(nrow, ncol))
}

convert_dt_to_sparse_matrix2 <- function(dt, i = 'i', j = 'j', x = 'value') {
  setorderv(dt, cols = c(j, i))
  mat <- data.table::dcast(dt, i ~ j, value.var = "value", fill = 0)
  rownames <- mat[[i]]
  mat <- as(mat[, -i], 'sparseMatrix')
  rownames(mat) <- rownames
  return(mat)
}

convert_matsbyname_to_dt <- function(x) {
  x <- as.matrix(x)
  x <- as.data.table(x, keep.rownames = 'i')
  melt(x, id.vars = 'i', variable.name = 'j', variable.factor = FALSE)

}



reshape_data <- function(x, use_names = TRUE) {
  if (inherits(x, "dgeMatrix")) {
    x <- as.matrix(x)
  }
  has_colnames <- !is.null(colnames(x))
  if (!has_colnames) {
    x <- as.data.table(x, keep.rownames = FALSE)
    x[, i := .I]
  } else {
    x <- as.data.table(x, keep.rownames = 'i')
  }

  x <- melt(x,
            id.vars = "i",
            variable.name = 'j')
  x <- x[abs(value) > 0]

  if (has_colnames) {
    x[, i := as.factor(i)]
  } else {
    x[, j := as.integer(j)]
  }

  return(x[])
}

cormat2sparse <- function(x, use_names = TRUE) {
  x[upper.tri(x, diag = TRUE)] <- NA
  return(reshape_data(x))
}



cordt <- function(x, method = 'pearson', use = 'everything',
                  nThreads = 6,
                  large = FALSE, use_names = TRUE) {
  if (inherits(x, 'data.frame')) x <- as.matrix(x)
  if (isTRUE(large))  cormat <- Rfast::cora(x, large = large)
  else cormat <- WGCNA::cor(x, method = method, use = use)

  return(cormat2sparse(cormat, use_names = use_names))
}

covdt <- function(x, method = 'pearson', use = 'everything',
                  nThreads = 6,
                  large = FALSE, use_names = TRUE) {
  if (inherits(x, 'data.frame')) x <- as.matrix(x)
  if (isTRUE(large))  cormat <- Rfast::cova(x, large = large)
  else cormat <- cov(x, method = method, use = use)

  return(cormat2sparse(cormat, use_names = use_names))
}

rlnorm_vec <- function(n, meanlog, sdlog) {
  return((mapply(rlnorm, meanlog = meanlog, sdlog = sdlog, n = n)))
}

rtruncnorm_vec <- function(n, a = rep(0, length(mean)), b = rep(Inf, length(mean)),
                           mean, sd) {
  return(mapply(rtruncnorm, n = n, a = a, mean = mean, sd = sd))
}


rlnorm2 <- function(n, mean = 1, sd = exp(1)) {

  # samples from lognormal with m1=1, s1=1:

  sig1 = sqrt(log(2))
  mu1  = - 1/2*log(2)

  y1 <- rlnorm(n, mu1, sig1)

  # then transforms the sample into a lognormal with m2=m, s2=2:

  sig2 = sqrt(log(1+sd^2/mean^2))
  mu2  = log(mean) - 1/2*sig2^2

  y2 <- exp(mu2)*(exp(-mu1)*y1)^(sig2/sig1)

  y2
}


rlnorm2_vec <- function(n, mean, sd) {
  return((mapply(rlnorm2, mean = mean, sd = sd, n = n)))
}

rnorm_vec <- function(n, mean, sd) {
  return((mapply(rnorm, mean = mean, sd = sd, n = n)))
}

#' like rgamma but parametrixed by mean and sd
#'
#' @param n
#' @param mean
#' @param sd
#'
#' @return
#' @export
#'
#' @examples
rgamma2 <- function(n, mean, sd) {
  shape <- (mean / sd) ^ 2
  rate <- mean / (sd) ^ 2
  rgamma(n, shape, rate)
}

rgamma2_vec <- function(n, mean, sd) {
  return(mapply(rgamma2, mean = mean, sd = sd, n = n))
}


#' from: https://github.com/r-forge/lcmix/blob/master/pkg/R/distributions.R
#'
#' @param n
#' @param corr correlation matrix
#' @param mean vector of means
#' @param sd vector of standard deviation
#'
#' @return
#' @export
#'
#' @examples
rmvgamma2 <- function(n, mean=10, sd=1, corr=diag(length(shape)))
{
  ## extract parameters, do sanity checks, deal with univariate case
  shape <-  (mean / sd) ^ 2
  rate <- mean / (sd) ^ 2
  rmvgamma(n, shape = shape, rate = rate, corr = corr)
}

# random generation function
rmvnorm <- function(n, mean=NULL, cov=NULL)
{
  ## munge parameters PRN and deal with the simplest univariate case

  if(is.null(mean))
    if(is.null(cov))
      return(rnorm(n))
  else
    mean = rep(0, nrow(cov))
  else if (is.null(cov))
    cov = diag(length(mean))

  ## gather statistics, do sanity checks

  D = length(mean)
  if (D != nrow(cov) || D != ncol(cov))
    stop("length of mean must equal nrow and ncol of cov")

  E = eigen(cov, symmetric=TRUE)
  if (any(E$val < 0))
    stop("Numerically negative definite covariance matrix")

  ## generate values and return

  mean.term = mean
  covariance.term = E$vec %*% (t(E$vec) * sqrt(E$val))
  independent.term = matrix(rnorm(n*D), nrow=D)

  drop(t(mean.term + (covariance.term %*% independent.term)))
}


# random generation function
rmvweisd <- function(n, shape=1, decay=1, corr=diag(length(shape)))
{
  ## extract parameters, do sanity checks, deal with univariate case

  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)

  Ds = length(shape)
  if(Ds > D)
    warning("'shape' longer than width of 'corr', truncating to fit")
  if(Ds != D)
    shape = rep(shape, length.out=D)

  Dd = length(decay)
  if(Dd > D)
    warning("'decay' longer than width of 'corr', truncating to fit")
  if(Dd != D)
    decay = rep(decay, length.out=D)

  if(D == 1) rweisd(n, shape, decay)

  ## generate standard multivariate normal matrix, convert to CDF

  Z = rmvnorm(n, cov=corr)
  cdf = pnorm(Z)

  ## convert to Weibull (WeiSD), return

  sapply(1:D, function(d) qweisd(cdf[,d], shape[d], decay[d]))
}


rmvweibull <- function(n, mean, sd, corr) {
  # see: https://stats.stackexchange.com/questions/159452/how-can-i-recreate-a-weibull-distribution-given-mean-and-standard-deviation-and
  #var <- sd^2
  shape <- (sd/mean)^(-1.086)
  scale <- mean / (gamma(1+(1/shape)))
  decay <- scale^(-shape)
  rmvweisd(n = n, shape = shape, decay = decay, corr = corr)
}

rweibull2 <- function(n, mean, sd) {
  # see: https://stats.stackexchange.com/questions/159452/how-can-i-recreate-a-weibull-distribution-given-mean-and-standard-deviation-and
  #var <- sd^2
  shape <- (sd/mean)^(-1.086)
  scale <- mean / (gamma(1+(1/shape)))
  decay <- scale^(-shape)
  rweibull(n, shape = shape, scale = scale)
}
unlist_sample <- function(x) {
  rbindlist(lapply(x, function(x) data.table('id_run' = 1:length(x),
                                             'value' = x)))
}


sparse.cor4 <- function(x){
  #https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r/9626089#9626089
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}
#' Compute correlation scores for columns of a sparse matrix
#'
#' @param m a sparse matrix, preferrably dgCMatrix
#'
#' @return a matrix of correlation values between each column of m
#'
sparse_cor <- function(m) {
  #library(Matrix)

  n_rows <- nrow(m)
  n_cols <- ncol(m)

  ii <- unique(m@i) + 1 # rows with a non-zero element

  Ex <- Matrix::colMeans(m)

  nozero <- as.vector(m[ii,]) - rep(Ex, each = length(ii))        # colmeans

  covmat <- (crossprod(matrix(nozero, ncol = n_cols)) +
               crossprod(t(Ex)) * (n_rows - length(ii))
  ) / (n_rows - 1)

  sdvec <- sqrt(diag(covmat))

  cormat <- covmat / crossprod(t(sdvec))

  return(cormat)
}

#' Compute covariance matrix for columns of a sparse matrix
#'
#' @param m a sparse matrix, preferrably dgCMatrix
#'
#' @return a matrix of correlation values between each column of m
#'
sparse_cov <- function(m) {
  #library(Matrix)

  n_rows <- nrow(m)
  n_cols <- ncol(m)

  ii <- unique(m@i) + 1 # rows with a non-zero element

  Ex <- Matrix::colMeans(m)

  nozero <- as.vector(m[ii,]) - rep(Ex, each = length(ii))        # colmeans

  covmat <- (crossprod(matrix(nozero, ncol = n_cols)) +
               crossprod(t(Ex)) * (n_rows - length(ii))
  ) / (n_rows - 1)


  return(covmat)
}
plot_dir_sample <- function(x, type = 'histogram') {
  x <- as.data.table(x)
  x <- suppressWarnings(melt(x))

  if (type == 'histogram') {
    ggplot(data = x, aes(x = value, col = variable, fill = variable)) +
      geom_histogram(alpha = 0.3, position = 'identity') +
      scale_fill_colorblind7() +
      scale_color_colorblind7()
  } else {
    stop('not implemented yet')
  }
}

summary2 <- function(x, na.rm = FALSE) {
  c(
    'mean' = mean(x, na.rm = na.rm),
    'sd' = sd(x, na.rm = na.rm)
  )
}


### Dirichlet functions ========================================================

#' Covariance matrix of Dirichlet random numbers
#'
#' @param alpha a vector of alpha parameter (shares * gamma)
#'
#' @return
#' @export
#'
#' @examples
cov_dirichlet <- function(alpha) {
  alpha0 <- sum(alpha)
  covmat <- matrix(0, nrow = length(alpha), ncol = length(alpha))

  diag(covmat) <- var_dirichlet(alpha)

  for (i in 1:length(alpha)) {
    for (j in 1:i) {
      if (i != j) covmat[i,j] <- -alpha[i] * alpha[j] / (alpha0^2 * (alpha0 + 1))
    }
  }
  covmat[upper.tri(covmat)] <- covmat[lower.tri(covmat)]
  return(covmat)
}

#' Correlation matrix of Dirichlet random nubmers
#'
#' @param alpha a vector of alpha parameter (shares * gamma)
#'
#' @return
#' @export
#'
#' @examples
cor_dirichlet <- function(alpha) {
  cov2cor(cov_dirichlet(alpha))
}


#' The variances of Dirichlet random numbers
#'
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
var_dirichlet <- function(alpha) {
  alpha0 <- sum(alpha)
  (alpha * (alpha0 - alpha)) / (alpha0^2 * (alpha0 + 1))
}

