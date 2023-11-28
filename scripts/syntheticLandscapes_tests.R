rm(list = ls())
source('./ecoFunctions.R')
source('./auxFunctions.R')

N <- 8 # number of species
comms <- makeComms(N)

# tune weight of higher-order interactions
scale_coeffs <- function(coeff_order, mean_scale = 1, sd_scale = 2) exp(-(coeff_order - mean_scale)^2/sd_scale^2)
scale_coeffs <- function(coeff_order, mean_scale = 1, sd_scale = 2) 0.2 + (1 - 0.2)*exp(-(coeff_order - mean_scale)^2/sd_scale^2)
scale_coeffs <- function(coeff_order, mean_scale = 1, sd_scale = 2) sapply(coeff_order,
                                                                           FUN = function(x) {
                                                                             
                                                                             if (x <= 1) return(1)
                                                                             else return(0.2 + (1 - 0.2)*2^(-1*(x - 1)/(sd_scale)))
                                                                             
                                                                           })

# interaction coefficients plots (wrapper function)
makeCoeffPlots <- function(sd_scale = 2) {
  
  mean_coeff <- rep(c(-1, 0, 1), 3)
  sd_coeff <- rep(c(0.1, 0.5, 0.9), each = 3)
  n_points <- 10000
  
  plot_this_all <- data.frame(case = character(0),
                              order = numeric(0),
                              y = numeric(0))
  for (i in 1:length(mean_coeff)) {
    
    plot_this <- data.frame(order = rep(1:N, each = n_points),
                            y = rnorm(n_points*N, mean = mean_coeff[i], sd = sd_coeff[i]))
    plot_this$y <- plot_this$y * scale_coeffs(plot_this$order, sd_scale = sd_scale)
    
    plot_this_all <- rbind(plot_this_all,
                           cbind(case = LETTERS[i],
                                 plot_this))
    
  }
  
  ylim <- 3
  plot_this_all <- plot_this_all[abs(plot_this_all$y) < ylim, ]
  
  myplot <-
    ggplot(plot_this_all,
           aes(x = order, y = y, group = order, fill = order)) +
      geom_abline(slope = 0,
                  intercept = 0,
                  color = 'gray',
                  linewidth = 0.5) +
      geom_violin(scale = 'width',
                  width = 0.5,
                  #fill = 'gray',
                  color = NA) +
      scale_y_continuous(name = expression(paste('Interaction coefficient, ', italic(delta), sep = '')),
                         limits = c(-ylim, ylim),
                         breaks = c(-1.5, 0, 1.5)) +
      scale_x_continuous(name = 'Order of interaction',
                         breaks = 1:N,
                         limits = c(1-0.5, N+0.5)) +
      facet_wrap(~case) +
      scale_fill_gradient(low = '#76d3d6', high = '#d32f37') +
      theme_bw() + theme(aspect.ratio = 0.4,
                         panel.grid = element_blank(),
                         axis.title = element_text(size = 16),
                         axis.text = element_text(size = 14),
                         panel.border = element_blank(),
                         legend.position = 'none',
                         strip.background = element_blank(),
                         strip.text = element_blank()) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
  return(myplot)
  
}

# coeff plots (small high-order epistasis)
#makeCoeffPlots(sd_scale = 2)

# coeff plots (small high-order epistasis)
#makeCoeffPlots(sd_scale = 4)

# coeff plots (large high-order epistasis)
#makeCoeffPlots(sd_scale = 100)


### GENERATE SYNTHETIC LANDSCAPES AND TEST PREDICTION METHOD


params_file <- list.files(path = './',
                          pattern = 'syntheticLandscapes_params_multi.txt',
                          full.names = T)

if (length(params_file)) {
  
  params <- read.table(params_file, sep = '\t', header = T)
  
} else {
  
  # params <- data.frame(mean_coeff = rep(seq(-1, 1, length.out = 7), 7),
  #                      sd_coeff = rep(seq(0.1, 0.9, length.out = 7), each = 7))
  # params <- do.call(rbind,
  #                   lapply(10^seq(log10(0.4), log10(15), length.out = 12),
  #                          FUN = function(x) cbind(params, sd_scale = x)))
  # params <- do.call(rbind, replicate(8, params, simplify = FALSE))
  params <- data.frame(mean_coeff = rep(seq(-1, 1, length.out = 7), 7),
                       sd_coeff = rep(seq(0.1, 0.9, length.out = 7), each = 7))
  params <- do.call(rbind,
                    lapply(10^seq(log10(0.4), log10(15), length.out = 12),
                           FUN = function(x) cbind(params, sd_scale = x)))
  params <- do.call(rbind, replicate(4, params, simplify = FALSE))
  
  params <- cbind(params,
                  rs = NA,
                  R2_stitching = NA,
                  R2.identity_stitching = NA,
                  relError.mean_topbot10_stitching = NA,
                  relError.sd_topbot10_stitching = NA,
                  R2_topbot10_stitching = NA,
                  R2.identity_topbot10_stitching = NA,
                  R2_reg1 = NA,
                  R2.identity_reg1 = NA,
                  relError.mean_topbot10_reg1 = NA,
                  relError.sd_topbot10_reg1 = NA,
                  R2_topbot10_reg1 = NA,
                  R2.identity_topbot10_reg1 = NA,
                  R2_reg2 = NA,
                  R2.identity_reg2 = NA,
                  relError.mean_topbot10_reg2 = NA,
                  relError.sd_topbot10_reg2 = NA,
                  R2_topbot10_reg2 = NA,
                  R2.identity_topbot10_reg2 = NA)
  
  write.table(as.data.frame(t(colnames(params))),
              file = './syntheticLandscapes_params_multi.txt',
              quote = F,
              row.names = F,
              col.names = F,
              sep = '\t')
  
  for (p in 1:nrow(params)) {
    
    print(paste('Synth. landscape ', p, ' of ', nrow(params), ' (', round(100*p/nrow(params)), '%)', sep = ''))
    
    # sample Fourier coefficients
    coeff <- setNames(unlist(lapply(0:N,
                                    FUN = function(i) apply(t(combn(N, i)),
                                                            FUN = function(x) rnorm(1,
                                                                                    mean = params$mean_coeff[p],
                                                                                    sd = params$sd_coeff[p]),
                                                            MARGIN = 1))),
                      unlist(lapply(0:N,
                                    FUN = function(i) apply(t(combn(N, i)),
                                                            FUN = function(x) paste(x, collapse = ''),
                                                            MARGIN = 1))))
    coeff <- coeff * scale_coeffs(nchar(names(coeff)), sd_scale = params$sd_scale[p])
    
    # make synthetic landscape
    synthLandscape <- cbind(comms,
                            fun = funFromCoeff(comms, coeff))
    #plotFitnessGraph(synthLandscape)
    #plotFEEs_clean(synthLandscape)
    
    # get ruggedness
    params$rs[p] <- get_rs(synthLandscape)
    
    # evaluate quality of predictions (stitching method) in synthetic landscape
    po <- evaluatePredictions_mult(synthLandscape)
    mylm <- lm(fun_true ~ fun_predicted, data = po)
    po_extremes <- po[po$fun_true > quantile(po$fun_true, probs = 0.9) | po$fun_true < quantile(po$fun_true, probs = 0.1), ]
    mylm_ext <- lm(fun_true ~ fun_predicted, data = po_extremes)
    
    params$R2_stitching[p] <- summary(mylm)$r.squared
    params$R2.identity_stitching[p] <- 1 - sum((po$fun_true - po$fun_predicted)^2)/sum((po$fun_true - mean(po$fun_true))^2)
    params$relError.mean_topbot10_stitching[p] <- mean(abs((po_extremes$fun_predicted - po_extremes$fun_true) / po_extremes$fun_true))
    params$relError.sd_topbot10_stitching[p] <- sd(abs((po_extremes$fun_predicted - po_extremes$fun_true) / po_extremes$fun_true))
    params$R2_topbot10_stitching[p] = summary(mylm_ext)$r.squared
    params$R2.identity_topbot10_stitching[p] = 1 - sum((po_extremes$fun_true - po_extremes$fun_predicted)^2)/sum((po_extremes$fun_true - mean(po_extremes$fun_true))^2)
    
    # evaluate quality of predictions by 1st and 2nd order regressions
    po <- get_all_loo_fits_mult(synthLandscape)
    
    # 1st order
    mylm <- lm(fun_true ~ fun_pred_1st, data = po)
    po_extremes <- po[po$fun_true > quantile(po$fun_true, probs = 0.9) | po$fun_true < quantile(po$fun_true, probs = 0.1), ]
    mylm_ext <- lm(fun_true ~ fun_pred_1st, data = po_extremes)
    
    params$R2_reg1[p] <- summary(mylm)$r.squared
    params$R2.identity_reg1[p] <- 1 - sum((po$fun_true - po$fun_pred_1st)^2)/sum((po$fun_true - mean(po$fun_true))^2)
    params$relError.mean_topbot10_reg1[p] <- mean(abs((po_extremes$fun_pred_1st - po_extremes$fun_true) / po_extremes$fun_true))
    params$relError.sd_topbot10_reg1[p] <- sd(abs((po_extremes$fun_pred_1st - po_extremes$fun_true) / po_extremes$fun_true))
    params$R2_topbot10_reg1[p] = summary(mylm_ext)$r.squared
    params$R2.identity_topbot10_reg1[p] = 1 - sum((po_extremes$fun_true - po_extremes$fun_pred_1st)^2)/sum((po_extremes$fun_true - mean(po_extremes$fun_true))^2)
    
    # 2nd order
    mylm <- lm(fun_true ~ fun_pred_2nd, data = po)
    po_extremes <- po[po$fun_true > quantile(po$fun_true, probs = 0.9) | po$fun_true < quantile(po$fun_true, probs = 0.1), ]
    mylm_ext <- lm(fun_true ~ fun_pred_2nd, data = po_extremes)
    
    params$R2_reg2[p] <- summary(mylm)$r.squared
    params$R2.identity_reg2[p] <- 1 - sum((po$fun_true - po$fun_pred_2nd)^2)/sum((po$fun_true - mean(po$fun_true))^2)
    params$relError.mean_topbot10_reg2[p] <- mean(abs((po_extremes$fun_pred_2nd - po_extremes$fun_true) / po_extremes$fun_true))
    params$relError.sd_topbot10_reg2[p] <- sd(abs((po_extremes$fun_pred_2nd - po_extremes$fun_true) / po_extremes$fun_true))
    params$R2_topbot10_reg2[p] = summary(mylm_ext)$r.squared
    params$R2.identity_topbot10_reg2[p] = 1 - sum((po_extremes$fun_true - po_extremes$fun_pred_2nd)^2)/sum((po_extremes$fun_true - mean(po_extremes$fun_true))^2)
    
    # save
    write.table(params[p, , drop = F],
                file = './syntheticLandscapes_params_multi.txt',
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t',
                append = T)
    
  }
  
}



### EMPIRICAL DATASETS

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

emp_ruggedness <- do.call(rbind,
                          lapply(1:length(data),
                                 FUN = function(i) {
                                   
                                   df <- data[[i]]
                                   colnames(df)[ncol(df)] <- 'fun'
                                   rs <- get_rs(df)
                                   
                                   return(data.frame(dataset = basename(files)[i],
                                                     rs = rs))
                                   
                                 }))

emp_ruggedness$dataset <- setNames(c('Bacterial starch hydrolysis',
                                'Bacterial butyrate secretion',
                                'Phytoplankton biomass',
                                'Above-ground plant biomass',
                                'Bacterial xylose oxidation',
                                'Bacterial pyoverdine secretion'),
                                basename(files)[1:6])[emp_ruggedness$dataset]
emp_ruggedness$dataset <- factor(emp_ruggedness$dataset, levels = c('Above-ground plant biomass',
                                                          'Phytoplankton biomass',
                                                          'Bacterial xylose oxidation',
                                                          'Bacterial starch hydrolysis',
                                                          'Bacterial butyrate secretion',
                                                          'Bacterial pyoverdine secretion'))





### MAKE PLOTS


# variance by order plots
focal_sd_scale <- unique(params$sd_scale) #[c(F, T, F, F, T, F, F, T, F, F, T, F)]
rel_var <- do.call(rbind,
                   lapply(focal_sd_scale,
                          FUN = function(x) {
                            
                            data.frame(sd_scale = x,
                                       coeff_order = 1:N,
                                       rel_var = scale_coeffs(1:N, sd_scale = x)^2 / sum(scale_coeffs(1:N, sd_scale = x)^2))
                            
                          }))

ggplot(rel_var, aes(x = coeff_order, y = 100*rel_var, fill = coeff_order)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~ sd_scale,
             nrow = 1,
             scales = 'free') +
  scale_fill_gradient(name = '',
                       low = '#d32f37',
                       high = '#76d3d6',
                       limits = c(0, N),
                       breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(name = 'Order of interaction',
                     limits = c(0.5, N + 0.5)) +
  scale_y_continuous(name = '% of functional\nvariance explained',
                     breaks = c(0, 50, 100),
                     limits = c(0, 100)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        panel.spacing = unit(1, "lines"),
        legend.position = 'none',
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/synthLandscapes/var_rel.pdf',
       device = 'pdf',
       dpi = 600,
       width = 150,
       height = 50,
       units = 'mm',
       limitsize = F)



#heatmaps
mycolors <- c('#99d7dc', '#176766', '#b33a3b')


tst <- aggregate(R2.identity_stitching ~ mean_coeff + sd_coeff + sd_scale,
                 data = params,
                 FUN = mean)

ggplot(tst,
       aes(x = mean_coeff, y = sd_coeff, fill = pmax(R2.identity_stitching, 0))) +
  geom_tile() +
  scale_y_continuous(name = expression(sigma[italic(beta)]),
                     expand = c(0, 0),
                     breaks = c(0.1, 0.5, 0.9)) +
  scale_x_continuous(name = expression(bar(beta)),
                     expand = c(0, 0),
                     breaks = c(-1, 0, 1)) +
  scale_fill_gradient2(name = expression(paste(italic(R)^2, ' predictions vs. observations')),
                       low = '#d32f37',
                       high = '#76d3d6',
                       limits = c(0, 1),
                       midpoint = 0.5,
                       breaks = pretty_breaks(n = 3)) +
  facet_wrap(~ sd_scale, nrow = 1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.spacing = unit(1, "lines"),
        legend.position = 'top') +
  guides(fill=guide_colorbar(ticks.colour = NA))

ggsave(filename = '../plots/synthLandscapes/R2_vs_structure_heatmaps.pdf',
       device = 'pdf',
       dpi = 600,
       width = 350,
       height = 125,
       units = 'mm',
       limitsize = F)



# R2 vs ruggedness
tst <- gather(params[, c('rs',
                         'R2.identity_stitching',
                         'R2.identity_reg1',
                         'R2.identity_reg2')],
              prediction_method, R2, 2:4, factor_key = FALSE)

# tst <- gather(params[, c('rs',
#                        'R2.identity_stitching',
#                        'R2.identity_reg1',
#                        'R2.identity_reg2')],
#               prediction_method, R2, 2:4, factor_key = FALSE)
tst$R2 <- pmax(0, tst$R2)
tst$rs <- log10(tst$rs)
tst$prediction_method <- setNames(c('1st order regression',
                                    '2nd order regression',
                                    'FEE concatenation'),
                                  unique(tst$prediction_method)[c(which(grepl('reg1', unique(tst$prediction_method))),
                                                                  which(grepl('reg2', unique(tst$prediction_method))),
                                                                  which(grepl('stitching', unique(tst$prediction_method))))])[tst$prediction_method]
tst$prediction_method <- factor(tst$prediction_method,
                                levels = c('1st order regression',
                                           '2nd order regression',
                                           'FEE concatenation'))




# fit logistic functions
mylns <- do.call(rbind,
                 lapply(unique(tst$prediction_method),
                        FUN = function(method) {
                          
                          lns_i <- nls(R2 ~ myLogistic(rs, A, B),
                                       tst[tst$prediction_method == method, ],
                                       start = list(A = 1, B = 0))
                          
                          return(cbind(prediction_method = method,
                                       as.data.frame(t(coef(lns_i)))))
                          
                        }))

ggplot(tst[tst$rs > -1 & tst$rs < 1, ],
       aes(x = rs, y = R2, color = prediction_method)) +
  geom_point(alpha = 0.05) +
  geom_vline(data = emp_ruggedness,
             aes(xintercept = log10(rs), color = dataset),
             linetype = 'dashed') +
  geom_function(fun = myLogistic, args = list(A = mylns$A[2], B = mylns$B[2]),
                color = mycolors[1],
                linewidth = 1) +
  geom_function(fun = myLogistic, args = list(A = mylns$A[3], B = mylns$B[3]),
                color = mycolors[2],
                linewidth = 1) +
  geom_function(fun = myLogistic, args = list(A = mylns$A[1], B = mylns$B[1]),
                color = mycolors[3],
                linewidth = 1) +
  scale_x_continuous(name = expression(paste('Ruggedness (', italic(r), '/', italic(s), ')')),
                     breaks = c(-1, 0, 1),
                     labels = c(expression(10^-1), expression(10^0), expression(10^1))) +
  scale_y_continuous(name = expression(paste(italic(R)^2, ' predictions vs. observations')),
                     breaks = c(0, 0.5, 1)) +
  scale_color_manual(name = '',
                     values = c(setNames(mycolors, levels(tst$prediction_method)),
                                setNames(c('#d6d62d',
                                           '#66b666',
                                           '#cb96c3',
                                           '#d72027',
                                           '#519ed7',
                                           'black'),
                                         levels(emp_ruggedness$dataset)))) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

ggsave(filename = '../plots/synthLandscapes/R2_vs_ruggedness.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 175,
       units = 'mm',
       limitsize = F)



tst2 <- tst[tst$rs > 0.5 & tst$rs < 0.75, ]

ggplot(tst2, aes(x = prediction_method, y = R2, fill = prediction_method)) +
  geom_violin(alpha = 1,
              width = 0.6,
              scale = 'width',
              color = NA) +
  scale_fill_manual(name = '',
                    values = mycolors) +
  scale_y_continuous(name = expression(italic(R)^2),
                     breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size = 16))

ggsave(filename = '../plots/synthLandscapes/R2_at_mediumrugedness.pdf',
       device = 'pdf',
       dpi = 600,
       width = 120,
       height = 100,
       units = 'mm',
       limitsize = F)
