rm(list = ls())

library(tidyverse)
library(glmnet)
library(cowplot)
source('ecoFunctions.R')

# load starch data
data <- read.csv('../data_sets/amyl_Sanchez-Gorostiaga2019.csv')
data <- aggregate(formula = function. ~ .,
                  data = data,
                  FUN = mean)

# species names
species <- setNames(c('B. cereus',
                      'B. megaterium',
                      'B. mojavensis',
                      'B. subtilis',
                      'B. thuringiensis',
                      'P. polymyxa'), colnames(data)[1:6])

# get functional effects
df <- makeGEdata(matrix2string(data))
df_bgs <- string2matrix(df[, c(1, 3)])
df <- cbind(as.data.frame(df_bgs[, 1:6]), df[, 2:4])








# predict functional effects using (global) FEEs
coefs.fees <- data.frame(intercept = numeric(0),
                         slope = numeric(0))
for (s in names(species)) {
  lmod <- lm(d_f ~ background_f,
             data = df[df$knock_in == s, ])
  coefs.fees <- rbind(coefs.fees,
                      data.frame(intercept = as.numeric(lmod$coefficients[1]),
                                 slope = as.numeric(lmod$coefficients[2])))
}
rownames(coefs.fees) <- names(species)

model.fees <- function(knock_in, background_composition, background_function) {
  coefs.fees[knock_in, 'intercept'] + coefs.fees[knock_in, 'slope']*background_function
}

df$predicted_df.fees <- NA
for (i in 1:nrow(df)) {
  df$predicted_df.fees[i] <- model.fees(df$knock_in[i], df[i, 1:6], df$background_f[i])
}










# predict functional effects using regression of interaction coefficients (1st order)

###
### THIS IS ALL ADAPTED FROM ABBY'S CODE (1ST ORDER REGRESSION)
###

N <- ncol(data) - 1 

# if (any(rowSums(df.s[,1:N]) == 0)){
#   df.s <- df.s[-which(rowSums(df.s[,1:N]) == 0),]
# }

communities <- sapply(1:nrow(data),
                      FUN = function(i) paste(names(species)[data[i, -ncol(data)] == 1], collapse = ','))
unique_communities <- unique(communities)
data <- data %>% mutate(community = communities)

y <- as.matrix(data$`function.`)   
f <- as.formula(y ~ .)
x <- model.matrix(f, data %>% dplyr::select(all_of(names(species))))

get_folds <- function(df, unique_communities, n_folds = 10){
  
  n_points <- length(unique_communities) 
  
  base_fold_size <- floor(n_points/n_folds)
  rem_fold_size <- floor(n_points %% n_folds)
  diff <- n_folds - rem_fold_size
  
  if (diff != n_folds){
    full_sample <- c(rep(seq(1:diff), base_fold_size), rep((diff + 1): n_folds, base_fold_size + 1))
  } else {
    full_sample <- rep(seq(1:n_folds), base_fold_size)
  }
  
  fold_sample <- sample(full_sample, size = length(full_sample), replace = FALSE)
  fold_data <- data.frame(community = unique_communities, fold_id = fold_sample)
  
  fold_id_data <- left_join(df, fold_data, by = 'community')
  fold_ids <- fold_id_data$fold_id
  
  return(fold_ids)
}

fold_ids <- get_folds(data, unique_communities, n_folds = 10)

cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids)

y_test <- as.matrix(data$`function.`)
f_test <- as.formula(y_test ~ .)
x_test <- model.matrix(f_test, data %>% dplyr::select(all_of(names(species))))
y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )

###
### END OF ABBY'S CODE
###

y_pred_cv <- setNames(as.numeric(y_pred_cv[, 1]), unique_communities)
names(y_pred_cv)[nchar(names(y_pred_cv)) == 0] <- 'empty'
names(y_pred_cv) <- orderName(names(y_pred_cv))

df$predicted_df.reg_1 <- NA
for (i in 1:nrow(df)) {
  
  bg <- colnames(df)[1:6]
  bg <- paste(bg[df[i, 1:6] == 1], collapse = ',')
  if (!nchar(bg)) bg <- 'empty'
  bg <- orderName(bg)
  
  bg_knockin <- setNames(df[i, 1:6], colnames(df)[1:6])
  bg_knockin[df$knock_in[i]] <- 1
  bg_knockin <- paste(colnames(bg_knockin)[bg_knockin == 1], collapse = ',')
  bg_knockin <- orderName(bg_knockin)
  
  df$predicted_df.reg_1[i] <- as.numeric(y_pred_cv[bg_knockin]) - as.numeric(y_pred_cv[bg])
  
}














# predict functional effects using regression of interaction coefficients (2nd order)

###
### THIS IS ALL ADAPTED FROM ABBY'S CODE (2ND ORDER REGRESSION)
###

N <- ncol(data) - 1 

# if (any(rowSums(df.s[,1:N]) == 0)){
#   df.s <- df.s[-which(rowSums(df.s[,1:N]) == 0),]
# }

communities <- sapply(1:nrow(data),
                      FUN = function(i) paste(names(species)[data[i, -ncol(data)] == 1], collapse = ','))
unique_communities <- unique(communities)
data <- data %>% mutate(community = communities)

y <- as.matrix(data$`function.`)   
f <- as.formula(y ~ .*.)
x <- model.matrix(f, data %>% dplyr::select(all_of(names(species))))

get_folds <- function(df, unique_communities, n_folds = 10){
  
  n_points <- length(unique_communities) 
  
  base_fold_size <- floor(n_points/n_folds)
  rem_fold_size <- floor(n_points %% n_folds)
  diff <- n_folds - rem_fold_size
  
  if (diff != n_folds){
    full_sample <- c(rep(seq(1:diff), base_fold_size), rep((diff + 1): n_folds, base_fold_size + 1))
  } else {
    full_sample <- rep(seq(1:n_folds), base_fold_size)
  }
  
  fold_sample <- sample(full_sample, size = length(full_sample), replace = FALSE)
  fold_data <- data.frame(community = unique_communities, fold_id = fold_sample)
  
  fold_id_data <- left_join(df, fold_data, by = 'community')
  fold_ids <- fold_id_data$fold_id
  
  return(fold_ids)
}

fold_ids <- get_folds(data, unique_communities, n_folds = 10)

cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids)

y_test <- as.matrix(data$`function.`)
f_test <- as.formula(y_test ~ .*.)
x_test <- model.matrix(f_test, data %>% dplyr::select(all_of(names(species))))
y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )

###
### END OF ABBY'S CODE
###

y_pred_cv <- setNames(as.numeric(y_pred_cv[, 1]), unique_communities)
names(y_pred_cv)[nchar(names(y_pred_cv)) == 0] <- 'empty'
names(y_pred_cv) <- orderName(names(y_pred_cv))

df$predicted_df.reg_2 <- NA
for (i in 1:nrow(df)) {
  
  bg <- colnames(df)[1:6]
  bg <- paste(bg[df[i, 1:6] == 1], collapse = ',')
  if (!nchar(bg)) bg <- 'empty'
  bg <- orderName(bg)
  
  bg_knockin <- setNames(df[i, 1:6], colnames(df)[1:6])
  bg_knockin[df$knock_in[i]] <- 1
  bg_knockin <- paste(colnames(bg_knockin)[bg_knockin == 1], collapse = ',')
  bg_knockin <- orderName(bg_knockin)
  
  df$predicted_df.reg_2[i] <- as.numeric(y_pred_cv[bg_knockin]) - as.numeric(y_pred_cv[bg])
  
}













# predict functional effects using FEEs (branched for polymyxa)
df$predicted_df.branched_fees <- NA
for (s in names(species)) {
  
  if (s == 'P') { # if the focal species is P. polymyxa...
    
    n <- df$`T` == 1 & df$knock_in == 'P' # these communities have B. thuringiensis in the background
    lmod <- lm(formula = d_f ~ background_f,
               data = df[n, ])
    df$predicted_df.branched_fees[n] <- as.numeric(lmod$coefficients[1] + lmod$coefficients[2]*df$background_f[n])
    
    n <- df$`T` == 0 & df$knock_in == 'P' # these communities DO NOT have B. thuringiensis in the background
    lmod <- lm(formula = d_f ~ background_f,
               data = df[n, ])
    df$predicted_df.branched_fees[n] <- as.numeric(lmod$coefficients[1] + lmod$coefficients[2]*df$background_f[n])
    
  } else { # if the focal species is NOT P. polymyxa...
    
    n <- df$`P` == 1 & df$knock_in == s # these communities have P. polymyxa in the background
    lmod <- lm(formula = d_f ~ background_f,
               data = df[n, ])
    df$predicted_df.branched_fees[n] <- as.numeric(lmod$coefficients[1] + lmod$coefficients[2]*df$background_f[n])
    
    n <- df$`P` == 0 & df$knock_in == s # these communities DO NOT have P. polymyxa in the background
    lmod <- lm(formula = d_f ~ background_f,
               data = df[n, ])
    df$predicted_df.branched_fees[n] <- as.numeric(lmod$coefficients[1] + lmod$coefficients[2]*df$background_f[n])
    
  }
  
}





























# plot
plot_this <- df[, c(7, 9:13)]
plot_this <- gather(plot_this, model, predicted_df, predicted_df.fees:predicted_df.branched_fees, factor_key=TRUE)
plot_this$model <- setNames(c('1st order regression',
                              '2nd order regression',
                              'FEEs',
                              'Branched FEEs'),
                            c('predicted_df.reg_1',
                              'predicted_df.reg_2',
                              'predicted_df.fees',
                              'predicted_df.branched_fees'))[as.character(plot_this$model)]
plot_this$model <- factor(plot_this$model, levels = c('1st order regression',
                                                      '2nd order regression',
                                                      'FEEs',
                                                      'Branched FEEs'))

plot_this <- plot_this[plot_this$model != 'FEEs', ]

xy.lims <- range(c(plot_this$d_f, plot_this$predicted_df))

ggplot(plot_this, aes(x = d_f, y = predicted_df, color = knock_in)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'gray') +
  geom_point(cex = 3) +
  facet_wrap(~ model,
             nrow = 1) +
  scale_x_continuous(name = expression(paste('Observed ', Delta, italic(F), sep = '')),
                     breaks = pretty_breaks(n = 3),
                     limits = xy.lims) +
  scale_y_continuous(name = expression(paste('Predicted ', Delta, italic(F), sep = '')),
                     breaks = pretty_breaks(n = 3),
                     limits = xy.lims) +
  scale_color_brewer(palette = 'Dark2',
                     name = 'Species',
                     labels = species) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12,
                                   face = 'italic'))

ggsave(filename = '../plots/frunfrys.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 150,
       units = 'mm',
       limitsize = F)

ggplot(plot_this, aes(x = d_f, y = predicted_df, color = knock_in)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'gray') +
  geom_point(cex = 3) +
  facet_grid(knock_in ~ model) +
  scale_x_continuous(name = expression(paste('Observed ', Delta, italic(F), sep = '')),
                     breaks = pretty_breaks(n = 3),
                     limits = xy.lims) +
  scale_y_continuous(name = expression(paste('Predicted ', Delta, italic(F), sep = '')),
                     breaks = pretty_breaks(n = 3),
                     limits = xy.lims) +
  scale_color_brewer(palette = 'Dark2',
                     name = 'Species',
                     labels = species) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12,
                                   face = 'italic'))

ggsave(filename = '../plots/frunfrys3.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 300,
       units = 'mm',
       limitsize = F)











# plot functional effects of B. cereus
plot_this <- df[df$knock_in == 'C', ]
plot_this$background_community <- matrix2string(plot_this[, c(1:6, 8)])$community
plot_this$background_community <- gsub(',', '', plot_this$background_community)

plot_this$background_community[plot_this$background_community == ''] <- 'Monoculture'
plot_this$background_community <- factor(plot_this$background_community,
                                         levels = c('Monoculture', plot_this$background_community[plot_this$background_community != 'Monoculture']))
plot_this$d_f[plot_this$background_community == 'Monoculture'] <- 2.2

ggplot(plot_this, aes(x = background_community, y = d_f, fill = background_community)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0,
             color = 'black') +
  scale_fill_manual(values = setNames(c('red', rep('gray', nrow(plot_this) - 1)),
                                      levels(plot_this$background_community))) +
  scale_x_discrete(name = 'Background consortium') +
  scale_y_continuous(name = expression(paste(Delta,italic(F), sep = ''))) +
  ggtitle('B. cereus') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        aspect.ratio = 0.3,
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14,
                                   angle = 45,
                                   hjust = 1),
        axis.title = element_text(size = 18),
        panel.background = element_blank(),
        plot.title = element_text(size = 16,
                                  face = 'italic'),
        legend.position = 'none')

ggsave(filename = '../plots/frunfrys2.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 150,
       units = 'mm',
       limitsize = F)














tst <- aggregate(formula = cbind(d_f, predicted_df.reg_1) ~ knock_in,
                 data = df,
                 FUN = mean)

xy.lims <- c(-5, 20)
 

ggplot(tst, aes(x = predicted_df.reg_1, y = d_f)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  scale_x_continuous(limits = xy.lims) +
  scale_y_continuous(limits = xy.lims) +
  theme_bw() +
  theme(aspect.ratio = 1)

