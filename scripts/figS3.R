rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(tidyverse)
library(cowplot)

# load data sets
files <- list.files('../pyoverdine_data', full.names = T)
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
    
    return(NA)
    
  } else {
    
    # make predictions and compare with observations
    eps <- inferAllResiduals(ge_data)
    predF <- predictF_fullClosure(oos$community,
                                  data,
                                  eps)
    po <- merge(predF, oos, by = 'community', suffixes = c('_pred', '_obs'))
    r_squared <- cor(po$fun_pred, po$fun_obs)^2
    
    return(r_squared)
    
  }

}

n_insample <- sort(rep(seq(60,200, by = 20), 200))
rsq <- rep(NA, length(n_insample))
for (i in 1:length(rsq)) {
  
  print(i)
  rsq[i] <- rsq_vs_nInSample(n_insample[i], data = data)
  
}

# plot
rsq <- data.frame(n_insample = n_insample,
                  rsq = rsq)
rsq <- do.call(data.frame,
               aggregate(formula = rsq ~ n_insample,
                         data = rsq,
                         FUN = function(x) c(mean = mean(x, na.rm = T), err_minus = quantile(x, 0.05, na.rm = T), err_plus = quantile(x, 0.95, na.rm = T))))
rsq$n_insample <- rsq$n_insample / 2^8

myplot <- 
  ggplot(rsq, aes(x = n_insample, y = rsq.mean, ymax = rsq.err_minus.5., ymin = rsq.err_plus.95.)) +
    geom_errorbar(width = 0) +
    geom_line() +
    geom_point(cex = 3) +
    scale_x_continuous(name = 'Fraction of communities in-sample',
                       breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(name = expression(paste('Predicted vs. observed ',italic(R^2))),
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c('0', '0.5', '1')) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

print(myplot)
ggsave(myplot,
       filename = '../plots/figS3.pdf',
       device = 'pdf',
       dpi = 600,
       width = 150,
       height = 90,
       units = 'mm',
       limitsize = F)











if(F) {

# wrapper function: leave part of the data as out_of_sample, fit FEEs with the remaining data, use method to predict function and compare to observed function
pred_vs_obs <- function(data, f_out_of_sample = 0.2, mode = 'full_closure') {
  
  # average replicates
  data <- matrix2string(data)
  data <- aggregate(formula = fun ~ .,
                    data = data,
                    FUN = mean)

  # leave 20% of data out-of-sample
  which_out_of_sample <- sample(1:nrow(data), size = round(f_out_of_sample*nrow(data)))
  out_of_sample <- data[which_out_of_sample, ]
  data <- data[-which_out_of_sample, ]
  
  # make sure the 'empty' community is left in sample
  empty_comm <- which(out_of_sample$community == '')
  if (length(empty_comm)) {
    data <- rbind(data, out_of_sample[empty_comm, ])
    out_of_sample <- out_of_sample[-empty_comm, ]
  }
  
  # fit FEEs with the remaining data
  ge_data <- makeGEdata(data)
  fits <- makeFEEs(ge_data)
  
  # make sure that all FEEs could be fitted (there is enough data in-sample for it), otherwise return NAs
  if (any(is.na(fits[, 'b']))) {
    
    return(list(df = data.frame(community = character(0),
                                fun = numeric(0),
                                predicted_f = numeric(0)),
                r_squared = numeric(0)))
    
  } else {
    
    if (mode == 'full_closure') {
    
      eps <- inferAllResiduals(ge_data)
      
      # use 'fullClosure' (global inference of epsilons) method to predict the function of the out of sample communities
      predicted_f <- predictF_fullClosure(out_of_sample$community,
                                          data,
                                          eps)
      
      pred_obs <- merge(out_of_sample, predicted_f, by = 'community', suffixes = c('_obs', '_pred'))
      
      r_squared <- cor(pred_obs$fun_obs, pred_obs$fun_pred)^2
      
      return(list(df = pred_obs,
                  r_squared = r_squared))
    
    } else if (mode == 'base') {
      
      # use base method (eps = 0) to predict the function of the out of sample communities
      predicted_f <- predictF_base(out_of_sample$community,
                                   data)
      
      pred_obs <- merge(out_of_sample, predicted_f, by = 'community', suffixes = c('_obs', '_pred'))
      
      r_squared <- cor(pred_obs$fun_obs, pred_obs$fun_pred)^2
      
      return(list(df = pred_obs,
                  r_squared = r_squared))
      
    }
    
  }
  
}

mode <- rep('full_closure', length(files))
mode[grepl('Clark', files)] <- 'base' # for the Clark et al. data, global inference is computationally out of reach

for (i in 1:5) {
  
  # for the phytoplankton biomass dataset (Ghedini et al., scale functions by 1e-4 for easier readability)
  if (i == 3) data[[i]][, ncol(data[[i]])] <- data[[i]][, ncol(data[[i]])]/1e4
  
  # get predicted vs observed plots and r_squared (repeat 500 times for every data set)
  po <- data.frame(run = numeric(0),
                   community = character(0),
                   fun_obs = numeric(0),
                   fun_pred = numeric(0))
  r_squared <- NULL
  
  for (n in 1:500) {
    print(c(i, n))
    po_i <- pred_vs_obs(data[[i]], mode = mode[[i]])
    
    po <- rbind(po, cbind(data.frame(run = rep(n, nrow(po_i$df))),
                          po_i$df))
    r_squared <- c(r_squared, po_i$r_squared)
  }
  
  po_i <- po_i$df
  range <- c(min(c(po_i$fun_obs, po_i$fun_pred)),
             max(c(po_i$fun_obs, po_i$fun_pred)))
  
  # plot
  plot1 <-
    ggplot(po_i, aes(x = fun_pred, y = fun_obs)) +
      geom_abline(slope = 1,
                  intercept = 0,
                  color = '#d1d3d4') +
      geom_point(shape = 1,
                 cex = 3) +
      scale_x_continuous(breaks = pretty_breaks(n = 3),
                         name = expression(paste('Predicted ', italic('F'), ' [a.u.]', sep = '')),
                         limits = range) +
      scale_y_continuous(breaks = pretty_breaks(n = 3),
                         name = expression(paste('Observed ', italic('F'), ' [a.u.]', sep = '')),
                         limits = range) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face = 'italic',
                                      size = 10),
            aspect.ratio = 0.6,
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = 'none') +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
  plot2 <-
    ggplot(data.frame(r = r_squared), aes(x = r)) +
      geom_histogram(aes(y = ..count../sum(..count..)),
                     bins = 15,
                     fill = 'black',
                     alpha = 0.75) +
      scale_x_continuous(limits = c(0, 1),
                         name = expression(italic(R^2)),
                         breaks = c(0, 0.5, 1),
                         labels = c('0', '0.5', '1')) +
      scale_y_continuous(name = '') +
      theme_bw() +
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face = 'italic',
                                      size = 10),
            aspect.ratio = 1,
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = 'none') +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
  myplot <- plot_grid(plot1, plot2,
                      align = 'v',
                      nrow = 1,
                      rel_heights = c(3/4, 1/4),
                      rel_widths = c(3/4, 1/4))
  
  print(myplot)
  ggsave(myplot,
         filename = paste('../plots/figS2_', i, '.pdf', sep = ''),
         device = 'pdf',
         dpi = 600,
         width = 200,
         height = 100,
         units = 'mm',
         limitsize = F)
  
}

}

