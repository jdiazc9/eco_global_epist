rm(list = ls())
source('./ecoFunctions.R')
source('auxFunctions.R')

# load data
files <- c(list.files('../data_sets', full.names = T),
           '../pyoverdine_data/training_set.csv', '../pyoverdine_data/test_set.csv')

data <- lapply(1:length(files),
               FUN = function(i) {
                 
                 df <- read.csv(files[i], header = T)
                 
                 if (i == 3) df$function. <- df$function. / 1e4
                 if (i %in% c(6, 7)) df <- cbind(df[, 1:8],
                                                 function. = rowMeans(df[, 9:ncol(df)]))
                 df <- aggregate(function. ~ .,
                                 data = df,
                                 FUN = mean)
                 
                 return(df)
                 
                })
data[[6]] <- rbind(data[[6]], data[[7]])
data <- data[1:6]

# wrapper function: fit self-consistent model
getSCfits <- function(df) {
  
  params0 <- c(rep(0, ncol(df) - 1), 1)
  gedf <- makeGEdata(matrix2string(df))
  
  # error function to minimize
  getSqErrors <- function(params) {
    
    names(params) <- c(unique(gedf$knock_in), 'r')
    
    model_df <- params[gedf$knock_in]*params['r'] + params[gedf$knock_in]*gedf$background_f
    absolute_error <- model_df - gedf$d_f
    sum_squared_error <- sum(absolute_error^2)
    
    return(sum_squared_error)
    
  }
  
  # minimize error function
  params_opt <- optim(params0, getSqErrors)
  
  # return fit parameters
  return(data.frame(knock_in = unique(gedf$knock_in),
                    slope = params_opt$par[1:(length(params_opt$par) - 1)],
                    intercept = params_opt$par[1:(length(params_opt$par) - 1)] * params_opt$par[length(params_opt$par)]))
  
}

# wrapper function: get predictions from self-consistent model via leave-one-out cross-validation
predictF_sc <- function(target_community, scfits, F0 = 0) {
  
  target_sp <- strsplit(target_community, split = ',')[[1]]
  f_target <- F0
  for (sp in target_sp) f_target <- f_target +
    scfits$intercept[scfits$knock_in == sp] +
    scfits$slope[scfits$knock_in == sp]*f_target
  
  return(f_target)
  
}

# get predictions vs observations
focal_dataset <- 6

df <- data[[focal_dataset]]
gedf <- makeGEdata(matrix2string(df))
scfits <- getSCfits(df)
po <- do.call(rbind,
              lapply(2:nrow(df),
                     FUN = function(i) {
                       
                       target_community <- orderName(paste(colnames(df[, 1:(ncol(df) - 1)])[df[i, 1:(ncol(df) - 1)] == 1],
                                                           collapse = ','))
                       
                       df_outofsample <- df[i, ]
                       df_insample <- df[-i, ]
                       
                       scfits <- getSCfits(df_insample)
                       fun_predicted <- predictF_sc(target_community, scfits)
                       
                       return(data.frame(community = target_community,
                                         fun_true = df_outofsample$function.,
                                         fun_predicted = fun_predicted))
                       
                     }))
mylm <- lm(fun_true ~ fun_predicted, data = po)

ggplot(gedf,
       aes(x = background_f, y = d_f)) +
  geom_abline(slope = 0, intercept = 0, color = 'gray') +
  geom_point() +
  geom_abline(data = scfits,
              aes(slope = slope, intercept = intercept, color = slope)) +
  scale_color_gradient2(low = 'firebrick1',
                        high = 'deepskyblue',
                        mid = 'black') +
  scale_x_continuous(name = 'Background function [a.u.]',
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = expression(paste(Delta*italic(F), ' [a.u.]')),
                     breaks = pretty_breaks(n = 2)) +
  facet_wrap(~ knock_in,
             nrow = 2) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        panel.border = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

# predicted vs observed
ggplot(po, aes(x = fun_predicted, y = fun_true)) +
  geom_blank(aes(x = fun_true, y = fun_predicted)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'gray') +
  geom_point() +
  scale_x_continuous(name = expression(paste(Predicted~italic(F), ' [a.u.]')),
                     breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(name = expression(paste(Empirical~italic(F), ' [a.u.]')),
                     breaks = pretty_breaks(n = 4)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = 'none')

ggsave(filename = '../plots/self_consistent_models/po_sc_pyoverdine.pdf',
       device = 'pdf',
       dpi = 600,
       width = 75,
       height = 75,
       units = 'mm',
       limitsize = F)

# predicted vs observed (inset)
ggplot(po, aes(x = fun_predicted, y = fun_true)) +
  #geom_blank(aes(x = fun_true, y = fun_predicted)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'gray') +
  geom_point() +
  scale_x_continuous(name = expression(paste(Predicted~italic(F), ' [a.u.]')),
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = expression(paste(Empirical~italic(F), ' [a.u.]')),
                     breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') +
  annotate("text",
           label = paste('R2 =', round(100 * summary(mylm)$r.squared)/100),
           x = min(po$fun_predicted) + 0.9*diff(range(po$fun_predicted)),
           y = min(po$fun_true) + 0.1*diff(range(po$fun_true)),
           size = 6, hjust = 1)

ggsave(filename = '../plots/self_consistent_models/po_sc_pyoverdine_inset.pdf',
       device = 'pdf',
       dpi = 600,
       width = 55,
       height = 55,
       units = 'mm',
       limitsize = F)

# compare that with the standard non-self-consistent model
po_standard <- evaluatePredictions(df)

mylm_standard <- lm(fun_true ~ fun_predicted, data = po_standard)

ggplot(po_standard, aes(x = fun_predicted, y = fun_true)) +
  geom_blank(aes(x = fun_true, y = fun_predicted)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'gray') +
  geom_point() +
  scale_x_continuous(name = expression(paste(Predicted~italic(F), ' [a.u.]')),
                     breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(name = expression(paste(Empirical~italic(F), ' [a.u.]')),
                     breaks = pretty_breaks(n = 4)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = 'none') +
  annotate("text",
           label = paste('R2 =', round(100 * summary(mylm_standard)$r.squared)/100),
           x = min(po_standard$fun_predicted) + 0.9*diff(range(po_standard$fun_predicted)),
           y = min(po_standard$fun_true) + 0.1*diff(range(po_standard$fun_true)),
           size = 6, hjust = 0.5)

ggsave(filename = '../plots/self_consistent_models/po_standard_pyoverdine.pdf',
       device = 'pdf',
       dpi = 600,
       width = 75,
       height = 75,
       units = 'mm',
       limitsize = F)
