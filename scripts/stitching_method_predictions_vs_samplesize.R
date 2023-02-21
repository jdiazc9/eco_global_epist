rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(tidyverse)
library(cowplot)

# load data sets
files <- list.files('../data_sets', full.names = T)
files <- files[c(1, 2, 5)] # for this analysis we exclude the 2 combinatorially smallest data sets
data <- lapply(files, FUN = function(file) read.csv(file))
data <- lapply(data, FUN = matrix2string)
data <- lapply(data,
               FUN = function(x) aggregate(formula = fun~community,
                                           data = x,
                                           FUN = mean))

# wrapper function: leave a fraction of communities out-of-sample, use method to predict from the in-sample ones, compare predictions and observations
rsq_vs_nInSample <- function(f, data = data, mode = 'full_closure') {
  
  # subsample data
  ge_data <- makeGEdata(data)
  ge_data.j <- lapply(unique(ge_data$knock_in),
                      FUN = function(x) {
                        
                        out <- ge_data[ge_data$knock_in == x, ]
                        out <- out[sample(1:nrow(out), size = max(ceiling(f*nrow(out)), 2), replace = F), ] # leave at least 2 points for each FEE (otherwise there is no fit possible)
                        
                        return(out)
                        
                      })
  ge_data.j <- do.call(rbind, ge_data.j)
  data.j <- data.frame(community = sapply(1:nrow(ge_data.j),
                                          FUN = function(i) {
                                            out <- c(strsplit(ge_data.j$background[i], split = ',')[[1]], ge_data.j$knock_in[i])
                                            out <- orderName(paste(out, collapse = ','))
                                            return(out)
                                          }),
                       fun = ge_data.j$background_f + ge_data.j$d_f)
  data.j <- rbind(data.j,
                  data.frame(community = ge_data.j$background,
                             fun = ge_data.j$background_f))
  data.j <- rbind(data.j,
                  data.frame(community = '',
                             fun = 0))
  data.j <- aggregate(formula = fun~community,
                      data = data.j,
                      FUN = mean)
  
  # how many communities are there in sample?
  f_data <- nrow(data.j)/nrow(data)
  
  # how many points were used per FEE? (average)
  n_fees <- sapply(unique(ge_data.j$knock_in), FUN = function(x) sum(ge_data.j$knock_in == x))
  n_fees <- mean(n_fees)
  
  # which communities are out of sample?
  which_oos <- data$community[!(data$community %in% data.j$community)]
  oos <- data[data$community %in% which_oos, ]
  
  # predict functions out of sample and compare with observations
  if (mode == 'full_closure') { # if possible, globally infer residuals
    
    eps <- inferAllResiduals(ge_data.j)
    predF <- predictF_fullClosure(oos$community,
                                  data.j,
                                  eps)
    
  } else if (mode == 'base') {# if global inference not possible, use base mode
    
    predF <- predictF_base(oos$community,
                           data.j)
    
  }
  
  if (nrow(predF) > 0) {
    
    po <- merge(predF, oos, by = 'community', suffixes = c('_pred', '_obs'))
    r_squared <- cor(po$fun_pred, po$fun_obs)^2
    
    best <- rbind(data.frame(community = predF$community,
                             fun = predF$fun,
                             source = 'pred'),
                  data.frame(community = data.j$community,
                             fun = data.j$fun,
                             source = 'obs'))
    best <- best[order(best$fun, decreasing = T), ]
    data.j <- rbind(data.j, oos)
    F_opt <- data.j$fun[data.j$community == best$community[1]]
    F_opt <- F_opt/max(data.j$fun)
    
  } else {
    
    r_squared <- NA
    F_opt <- NA
    
  }

  
  return(c(f_data = f_data, n_fees = n_fees, r_squared = r_squared, F_opt = F_opt))
  
}

# apply wrapper function to all data sets
f <- sort(runif(n = 100, min = 0.1, max = 0.5)) # fraction of points kept in-sample to fit FEEs
mode <- rep('full_closure', length(data))
mode[grepl('Clark', files)] <- 'base' # we globaly infer residuals unless computationally out of reach

myruns <- data.frame(dataset = numeric(0),
                     f_sample = numeric(0),
                     n_fees = numeric(0),
                     r_squared = numeric(0),
                     F_opt = numeric(0))

for (i in 1:length(data)) {
  for (j in 1:length(f)) {

    print(paste('Dataset', i, '/ iteration', j))
    run <- rsq_vs_nInSample(f[j], data = data[[i]], mode = mode[i])
    
    myruns <- rbind(myruns,
                    data.frame(dataset = i,
                               f_sample = run['f_data'],
                               n_fees = run['n_fees'],
                               r_squared = run['r_squared'],
                               F_opt = run['F_opt']))
    
  }
}

# plot
plot_this <- gather(myruns, metric, value, r_squared:F_opt)
plot_this$metric <- factor(plot_this$metric, levels = c('r_squared', 'F_opt'))

ggplot(plot_this, aes(x = n_fees, y = value, color = metric)) +
  geom_smooth(method = "glm",
              formula = y ~ x,
              method.args = list(family = "binomial"), 
              se = TRUE,
              size = 0.5) +
    geom_point() +
    facet_wrap(metric ~ dataset,
               scales = 'free') +
    scale_x_continuous(name = 'Avg. number of\ndata points per fitted FEE',
                       breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(name = expression(paste(italic(R)^2, ' predicted vs. observed ', italic(F))),
                       limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c('0', '0.5', '1')) +
    scale_color_manual(values = c('black', 'red')) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 18),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.text.y.right = element_text(size = 16,
                                           color = 'red'),
          axis.title = element_text(size = 18),
          axis.title.y.right = element_text(size = 18,
                                            color = 'red'),
          axis.ticks.y.right = element_line(color = 'red'),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/stitching_rsquared_Fopt_vs_samplesize.pdf',
       device = 'pdf',
       dpi = 600,
       width = 250,
       height = 180,
       units = 'mm',
       limitsize = F)

