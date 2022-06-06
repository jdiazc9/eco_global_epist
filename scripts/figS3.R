rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(tidyverse)
library(cowplot)

# load data sets
files <- list.files('../pyoverdine_data', full.names = T, recursive = F, pattern = '.csv')
data <- lapply(files, FUN = function(file) read.csv(file))

data <- lapply(data,
               FUN = function(data) cbind(data[, 1:8], fun = rowMeans(data[, 9:ncol(data)])))
data <- do.call(rbind, data)
data <- matrix2string(data)

# wrapper function: leave a fraction of communities out-of-sample, use method to predict from the in-sample ones, compare predictions and observations
rsq_vs_nInSample <- function(n_in_sample, data = data) {
  
  which_oos <- sample(1:nrow(data), size = nrow(data) - n_in_sample, replace = F)
  oos <- data[which_oos, ]
  data <- data[-which_oos, ]
  
  # make sure the 'empty' community is left in sample
  empty_comm <- which(oos$community == '')
  if (length(empty_comm)) {
    data <- rbind(data, oos[empty_comm, ])
    oos <- oos[-empty_comm, ]
  }
  
  # make FEEs
  ge_data <- makeGEdata(data)
  fits <- makeFEEs(ge_data)
  
  # make sure all FEEs were fitted (enough points in-sample), otherwise return NA
  if (any(is.na(fits[, 'b']))) {
    
    return(c(r_squared = NA, F_opt = NA))
    
  } else {
    
    # make predictions and compare with observations
    eps <- inferAllResiduals(ge_data)
    predF <- predictF_fullClosure(oos$community,
                                  data,
                                  eps)
    po <- merge(predF, oos, by = 'community', suffixes = c('_pred', '_obs'))
    r_squared <- cor(po$fun_pred, po$fun_obs)^2
    
    # predicted best community
    best <- rbind(data.frame(community = predF$community,
                             fun = predF$fun,
                             source = 'pred'),
                  data.frame(community = data$community,
                             fun = data$fun,
                             source = 'obs'))
    best <- best[order(best$fun, decreasing = T), ]
    
    data <- rbind(data, oos)
    F_opt <- data$fun[data$community == best$community[1]]
    F_opt <- F_opt/max(data$fun)
    
    return(c(r_squared = r_squared, F_opt = F_opt))
    
  }

}

n_insample <- sort(rep(seq(60, 200, by = 20), 500))
rsq <- rep(NA, length(n_insample))
F_opt <- rep(NA, length(n_insample))
for (i in 1:length(rsq)) {
  
  print(i)
  run <- rsq_vs_nInSample(n_insample[i], data = data)
  rsq[i] <- run['r_squared']
  F_opt[i] <- run['F_opt']
  
}

# plot
rsq <- data.frame(n_insample = n_insample,
                  rsq = rsq,
                  F_opt = F_opt)
rsq <- do.call(data.frame,
               aggregate(formula = cbind(rsq, F_opt) ~ n_insample,
                         data = rsq,
                         FUN = function(x) c(mean = mean(x, na.rm = T), err_minus = quantile(x, 0.05, na.rm = T), err_plus = quantile(x, 0.95, na.rm = T))))
rsq$n_insample <- rsq$n_insample / 2^8

rsq <- rbind(data.frame(n_insample = rsq$n_insample,
                        metric = 'R_squared',
                        mean = rsq$rsq.mean,
                        err_minus = rsq$rsq.err_minus.5.,
                        err_plus = rsq$rsq.err_plus.95.),
             data.frame(n_insample = rsq$n_insample,
                        metric = 'F_opt',
                        mean = rsq$F_opt.mean,
                        err_minus = rsq$F_opt.err_minus.5.,
                        err_plus = rsq$F_opt.err_plus.95.))
rsq$n_insample[rsq$metric == 'R_squared'] <- rsq$n_insample[rsq$metric == 'R_squared'] + 0.005
rsq$n_insample[rsq$metric == 'F_opt'] <- rsq$n_insample[rsq$metric == 'F_opt'] - 0.005

myplot <- 
  ggplot(rsq, aes(x = n_insample, y = mean, ymin = err_minus, ymax = err_plus, color = metric)) +
    geom_errorbar(width = 0) +
    geom_line() +
    geom_point(cex = 3) +
    scale_x_continuous(name = 'Fraction of communities in sample',
                       breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(name = expression(paste('Predicted vs. observed ',italic(R^2))),
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c('0', '0.5', '1'),
                       sec.axis = sec_axis(trans=~.*1,
                                           name = expression(paste(italic(F), ' ( ', bold(s)[opt], ' ) / ', italic(F)[max], sep = '')),
                                           breaks = c(0, 0.5, 1),
                                           labels = c('0', '0.5', '1'))) +
    scale_color_manual(values = c('red', 'black')) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.text.y.right = element_text(size = 16,
                                           color = 'red'),
          axis.title = element_text(size = 18),
          axis.title.y.right = element_text(size = 18,
                                            color = 'red'),
          axis.ticks.y.right = element_line(color = 'red'),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
    annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5, color = 'red')

print(myplot)
ggsave(myplot,
       filename = '../plots/figS3.pdf',
       device = 'pdf',
       dpi = 600,
       width = 150,
       height = 90,
       units = 'mm',
       limitsize = F)

# plot insets
set.seed(0)

data_i <- data[sample(1:nrow(data), size = min(n_insample), replace = F), ]
ge_data <- makeGEdata(data_i)

myplot <-
  ggplot(ge_data[ge_data$knock_in == 'sp_1', ], aes(x = background_f, y = d_f)) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_point(color = 'black',
               cex = 3,
               shape = 16) +
    geom_smooth(method = 'lm',
                formula = y~x,
                se = FALSE,
                fullrange = TRUE,
                color = 'firebrick1') +
    scale_x_continuous(name = expression(paste(italic(F), '(background) [uM]')),
                       breaks = pretty_breaks(n = 2),
                       limits = c(0, 80)) +
    scale_y_continuous(name = 'dF [uM]',
                       breaks = pretty_breaks(n = 2),
                       limits = c(-35, 20)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
print(myplot)
ggsave(myplot,
       filename = '../plots/figS3_inset1.pdf',
       device = 'pdf',
       dpi = 600,
       width = 80,
       height = 50,
       units = 'mm',
       limitsize = F)

data_i <- data[sample(1:nrow(data), size = max(n_insample), replace = F), ]
ge_data <- makeGEdata(data_i)

myplot <- 
  ggplot(ge_data[ge_data$knock_in == 'sp_1', ], aes(x = background_f, y = d_f)) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_point(color = 'black',
               cex = 3,
               shape = 16) +
    geom_smooth(method = 'lm',
                formula = y~x,
                se = FALSE,
                fullrange = TRUE,
                color = 'firebrick1') +
    scale_x_continuous(name = expression(paste(italic(F), '(background) [uM]')),
                       breaks = pretty_breaks(n = 2),
                       limits = c(0, 80)) +
    scale_y_continuous(name = 'dF [uM]',
                       breaks = pretty_breaks(n = 2),
                       limits = c(-35, 20)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
print(myplot)
ggsave(myplot,
       filename = '../plots/figS3_inset2.pdf',
       device = 'pdf',
       dpi = 600,
       width = 80,
       height = 50,
       units = 'mm',
       limitsize = F)
