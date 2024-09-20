library(disaggR)
library(data.table)
library(tidyverse)
library(corrplot)
library(patchwork)
library(GGally)
source('functions.R')

my_fn <- function(data, mapping, method="p", use="pairwise", ...){

  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  # calculate correlation
  corr <- cor(x, y, method=method, use=use)

  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue,
  # zero to white, and one to red
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("red", "white", "blue"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]

  ggally_cor(data = data, mapping = mapping, ...) +
    theme(panel.background = element_rect(fill=fill, colour=NA),
          panel.grid.major = element_blank())
}


N <- 1E3

sample_shares <- list()
sample_shares[[1]] <- rshares(N, c(0.1, 0.3, 0.6))
sample_shares[[2]] <- rshares(N, c(0.1, 0.3, 0.6), sds = 0.18 * c(0.1, 0.3, 0.6))
sample_shares[[3]] <- rshares(N, c(0.1, 0.3, 0.6), sds = 0.01 * c(0.1, 0.3, 0.6))

sample_agg <- list()
sample_agg[[1]] <- ragg(N, mean= 100, sd = 1)
sample_agg[[2]] <- ragg(N, mean= 100, sd = 10)
sample_agg[[3]] <- ragg(N, mean= 100, sd = 10)

sample_disagg <- Map(`*`, x = sample_shares, y = sample_agg)

cors <- lapply(sample_disagg, cor)
cors

# for (i in 1:length(cors)) {
#   pdf(height=7, width=7, file=paste0("figures_overleaf/figures/correlation_generic",
#                                      i, ".pdf"))
#   corrplot(cors[[i]], diag = FALSE, type = 'lower', tl.col = 'black', tl.cex = 2)
#   dev.off()
#
# }


plist <- lapply(sample_disagg, function(x) {
  #sample_disagg[[1]] %>% as.data.table
  ggpairs(as.data.table(x) %>% setnames(as.character(1:3)),
          upper = list(continuous = function(...) my_fn(..., col = 'black',
                                                        title = 'Corr', digits = 2,
                                                        stars = FALSE)),
          #lower = list(continuous = wrap('points', )),
          lower = list(continuous = function(...) ggally_points(..., alpha = 0.4,
                                                                col = 'grey30',
                                                                shape = 16)+
                         theme_bw()),
          diag = list(continuous = function(...) ggally_barDiag(..., col = 'grey30',
                                                                fill = 'grey70')+
                        theme_bw())
  )

})

plist2 <- list()
for (i in 1:length(plist)) {
  filename <- paste0("figures_overleaf/figures/figureS1_",
                     i, ".pdf")
  ggsave(plot = plist[[i]],
         filename = filename,
         width = 4, height = 4)

}



