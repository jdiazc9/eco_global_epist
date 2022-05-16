rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(tidyverse)
library(cowplot)

# set seed for reproducibility
set.seed(0)

# load data sets
data <- lapply(list.files('../data_sets/', full.names = T), FUN = function(file) read.csv(file))

# wrapper function: leave part of the data as out_of_sample, fit FEEs with the remaining data, use method to predict function and compare to observed function
pred_vs_obs <- function(data, f_out_of_sample = 0.2) {
  
  # leave 20% of data out-of-sample (make sure that the 'empty' community is not one of them)
  which_out_of_sample <- sample(1:nrow(data), size = round(0.2*nrow(data)))
  out_of_sample <- data[which_out_of_sample, ]
  data <- data[-which_out_of_sample, ]
  empty_comm <- which(sapply(1:nrow(out_of_sample),
                             FUN = function(i) all(out_of_sample[i, 1:(ncol(out_of_sample) - 1)] == 0)))
  if (length(empty_comm) >  0) {
    data <- rbind(data, out_of_sample[empty_comm, ])
    out_of_sample <- out_of_sample[-empty_comm, ]
  }
  
  # fit FEEs with the remaining data
  ge_data <- makeGEdata(matrix2string(data), exclude.single.mut = F)
  fits <- makeGEfits(ge_data)
  eps <- inferEps(ge_data)
  
  # use method to predict the function of the out of sample communities
  out_of_sample <- matrix2string(out_of_sample)
  predicted_f <- predictF(out_of_sample$community,
                          data = matrix2string(data),
                          coeff = fits,
                          eps = eps)
  predicted_f <- data.frame(community = names(predicted_f),
                            predicted_f = as.numeric(predicted_f))
  
  pred_obs <- merge(out_of_sample, predicted_f)
  
  r_squared <- cor(pred_obs$fun, pred_obs$predicted_f)^2
  
  return(list(df = pred_obs,
              r_squared = r_squared))
  
}

i <- 3 #for (i in 1:5) {

  # if there are multiple measurements of a same community, average them
  colnames(data[[i]])[ncol(data[[i]])] <- 'fun'
  data[[i]] <- aggregate(formula = fun ~ .,
                         data = data[[i]],
                         FUN = mean)
  
  # for the phytoplankton biomass dataset (Ghedini et al., scale functions by 1e-4 for easier readability)
  if (i == 3) data[[i]][, ncol(data[[i]])] <- data[[i]][, ncol(data[[i]])]/1e4
  
  # get predicted vs observed plots and r_squared (repeat 50 times for every data set)
  po <- data.frame(community = character(0),
                   fun = numeric(0),
                   predicted_f = numeric(0))
  r_squared <- NULL
  
  for (n in 1:50) {
    print(c(i, n))
    po_i <- pred_vs_obs(data[[i]])
    
    po <- rbind(po, po_i$df)
    r_squared <- c(r_squared, po_i$r_squared)
  }
  
  po_i <- po_i$df
  
  range <- c(min(c(po_i$predicted_f, po_i$fun)),
             max(c(po_i$predicted_f, po_i$fun)))
  
  # plot
  plot1 <-
    ggplot(po_i, aes(x = predicted_f, y = fun)) +
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
         filename = '../plots/figS2.pdf',
         device = 'pdf',
         dpi = 600,
         width = 200,
         height = 100,
         units = 'mm',
         limitsize = F)
  
#}


