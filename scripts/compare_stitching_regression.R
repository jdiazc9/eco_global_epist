rm(list = ls())

library(tidyverse)
library(glmnet)
library(cowplot)

source('ecoFunctions.R')

##########################
####   functions    ######
###########################

#metric 
get_rmse <- function(observed, predicted){
  return(sqrt(sum((observed - predicted)^2)/length(observed)))
}

#helper function for cross-validated regression
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

####################################################
# get leave-one-out fit for a specified dataset #### 
#####################################################

get_all_loo_fits <- function(df, v = FALSE){
  
  N <- ncol(df) -1 
  
  all_species <- colnames(df)[1:N]
  communities <- sapply(1:nrow(df),
                        FUN = function(i) paste(all_species[df[i, -ncol(df)] == 1], collapse = ','))
  unique_communities <- unique(communities)
  df <- df %>% mutate(community = communities) 
  
  #add mean fitness to df 
  df <- df %>% group_by(community) %>% mutate(mean_fitness = mean(`function`)) %>% ungroup() 
  
  mean_df <- df %>% distinct(community, .keep_all = TRUE) %>% dplyr::select(-c(fitness))
  
  ####################
  #### stitching #####
  ####################
  
  loo_res <- tibble()
  for (exp in unique_communities){
    exp_inds <- which(mean_df$community == exp)
    
    data <- mean_df[-exp_inds,] %>% dplyr::select(-community) 
    if (! any(rowSums(mean_df[,1:N]) == 0)){
      data <- rbind(data, rep(0, (N + 1)))
    }
    
    obsF <- mean_df[exp_inds[1],] %>% dplyr::select(-community) 
    target <- matrix2string(obsF)
    
    data <- matrix2string(data)
    ge_data <- makeGEdata(data)
    
    eps <- try(inferAllResiduals(ge_data))
    
    predF <- predictF_fullClosure(target$community, data, eps)
    po <- merge(target, predF, by = 'community', suffixes = c('_obs', '_pred'))
    colnames(po) <- c('community', 'observed', 'predicted')
    
    tmp <- list(obs = po$observed, pred = po$predicted)
    loo_res <- rbind(loo_res, tmp)
  }
  
  #################################
  # first order linear regression #
  #################################
  loo_res_first_order <- tibble()
  for (exp in unique_communities){
    exp_inds <- which(df$community == exp)
    #set up regression
    y <- as.matrix(df[-exp_inds,]$fitness)   
    f <- as.formula(y ~ .)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$fitness)
    f_test <- as.formula(y_test ~ .)
    x_test <- model.matrix(f_test, df[exp_inds[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
    
    tmp <- list(obs = df[exp_inds[1],]$mean_fitness, pred = y_pred_cv)
    loo_res_first_order <- rbind(loo_res_first_order, tmp)
  }
  
  
  #################################
  # second order linear regression #
  #################################
  count <- 1
  loo_res_second_order <- tibble()
  for (exp in unique_communities){
    if (v){
      print(paste0(count , ' out of ', length(unique_communities)))
      count <- count + 1
    }
    
    exp_inds <- which(df$community == exp)
    #set up regression
    y <- as.matrix(df[-exp_inds,]$fitness)   
    f <- as.formula(y ~ .*.)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$fitness)
    f_test <- as.formula(y_test ~ .*.)
    x_test <- model.matrix(f_test, df[exp_inds[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
    
    tmp <- list(obs = df[exp_inds[1],]$mean_fitness, pred = y_pred_cv)
    loo_res_second_order <- rbind(loo_res_second_order, tmp)
  }
  
  r2 <- cor(loo_res$obs, loo_res$pred)^2
  r2_first_order <- cor(loo_res_first_order$obs, loo_res_first_order$pred)^2
  r2_second_order <- cor(loo_res_second_order$obs, loo_res_second_order$pred)^2
  
  return(list(r2 = r2, r2_first_order = r2_first_order, r2_second_order = r2_second_order, 
              loo_res = loo_res, loo_res_first_order = loo_res_first_order,
              loo_res_second_order = loo_res_second_order))
}

#######################################
###     get model fits with         ###
###    randomized % out of sample   ###
#######################################

fit_all_models_oof <- function(df, N, reps, v = FALSE){
  res <- tibble()
  
  all_species <- colnames(df)[1:N]
  communities <- sapply(1:nrow(df),
                        FUN = function(i) paste(all_species[df[i, -ncol(df)] == 1], collapse = ','))
  unique_communities <- unique(communities)
  df <- df %>% mutate(community = communities) 
  
  #add mean fitness to df 
  df <- df %>% group_by(community) %>% mutate(mean_fitness = mean(fitness)) %>% ungroup() 
  
  mean_df <- df %>% distinct(community, .keep_all = TRUE) %>% dplyr::select(-c(fitness))
  
  infit_range <- c(.6, .7, .8, .9)
  for (infit_pct in infit_range){
    for (i in 1:reps){
      if (v == TRUE){print(paste0('infit pct is: ', infit_pct, 'rep: ', i))} 
      
      #get random training sample of communities  
      train_ids <- sample(1:length(unique_communities), (infit_pct) * length(unique_communities))
      train <- which(df$community %in% unique_communities[train_ids])
      
      train_stitching <- which(mean_df$community %in% unique_communities[train_ids])
      data <- mean_df[train_stitching,] %>% dplyr::select(-community)
      
      #need to add 0 community
      data <- rbind(data, rep(0, (N + 1)))
      #process data for stitching procedure
      obsF <- mean_df[-train_stitching,] %>% dplyr::select(-community)
      target <- matrix2string(obsF)
      
      data <- matrix2string(data)
      ge_data <- makeGEdata(data)
      
      eps <- try(inferAllResiduals(ge_data))
      
      #get stitching predictions
      predF <- predictF_fullClosure(target$community, data, eps)
      po <- merge(target, predF, by = 'community', suffixes = c('_obs', '_pred'))
      colnames(po) <- c('community', 'observed', 'predicted')
      
      #get stitching metrics 
      r2 <- cor(po$observed, po$predicted)^2
      rmse_stitching <- get_rmse(po$observed, po$predicted)
      
      #################################
      # first order linear regression #
      #################################
      
      #set up regression
      y <- as.matrix(df[train,]$fitness)   
      f <- as.formula(y ~ .)
      x <- model.matrix(f, df[train,] %>% dplyr::select(all_of(all_species)))
      
      unique_communities_train <- unique_communities[train_ids]
      fold_ids <- get_folds(df[train,], unique_communities_train, n_folds = 10)
      
      cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
      
      #get out of fit cv 
      y_test <- as.matrix(df[-train,]$fitness)
      f_test <- as.formula(y_test ~ .)
      x_test <- model.matrix(f_test, df[-train,] %>% dplyr::select(all_of(all_species)))
      y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
      
      #bind predictions
      all_pred <- df[-train,] %>% mutate(predicted = y_pred_cv[,1]) %>% 
        distinct(community, .keep_all = TRUE)
      
      #get metrics cv
      r2_lin <- cor(all_pred$mean_fitness, all_pred$predicted)^2 
      rmse_lin <- get_rmse(all_pred$mean_fitness, all_pred$predicted)
      
      
      #################################
      # second order linear regression #
      #################################
      y <- as.matrix(df[train,]$fitness)
      f <- as.formula(y ~ .*.)
      x <- model.matrix(f, df[train,] %>% dplyr::select(all_of(all_species)))
      
      unique_communities_train <- unique_communities[train_ids]
      fold_ids <- get_folds(df[train,], unique_communities_train, n_folds = 10)
      
      cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
      
      #get out of fit cv 
      y_test <- as.matrix(df[-train,]$fitness)
      f_test <- as.formula(y_test ~ .*.)
      x_test <- model.matrix(f_test, df[-train, 1:N])
      y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
      
      #bind predictions
      all_pred <- df[-train,] %>% mutate(predicted = y_pred_cv[,1]) %>% 
        distinct(community, .keep_all = TRUE)
      
      #get metrics cv
      r2_cv_pairs <- cor(all_pred$mean_fitness, all_pred$predicted)^2 
      rmse_cv_pairs <- get_rmse(all_pred$mean_fitness, all_pred$predicted)
      
      plot(po$observed, po$predicted, main = 'stitching pred')
      abline(0,1)
      
      plot(all_pred$mean_fitness, all_pred$predicted, main = 'cv pred')
      abline(0,1)
      
      #bind this all together 
      tmp <- list(r2 = r2, r2_lin = r2_lin, r2_cv_pairs = r2_cv_pairs, 
                  rmse_stitching = rmse_stitching, rmse_lin = rmse_lin,
                  rmse_cv_pairs = rmse_cv_pairs, 
                  infit_pct = infit_pct)
      res <- rbind(res, tmp)
    }
  }
  return(res)
}

#################################
### get leave-one-out fits for ##
########  all datasets  #########
#################################

#list all dataframe paths
path <- '../data_sets/'
all_datasets <- list.files(path = path)

# want GE_biomass, Kuebbing natives, langenheder, pyoverdine, starch
all_datasets <- all_datasets[c(2,3, 6,7,8,9)]

p_list <- list()
res <- tibble() 
for (dataset_id in 1:length(all_datasets)){
  df <- read_csv(paste0(path, all_datasets[dataset_id]), show_col_types = F)
  
  name <- str_remove(all_datasets[dataset_id], '.csv')
  
  if (name == 'starch_data'){
    df <- df[-which(rowSums(df[,1:6]) == 0),]
  }
  
  loo_fits <- get_all_loo_fits(df) 
  
  tmp <- list(name = name, r2 = loo_fits$r2, r2_first_order = loo_fits$r2_first_order,
              r2_second_order = loo_fits$r2_second_order)
  res <- rbind(res, tmp)
  
  max_scale <- max(max(loo_fits$loo_res$obs, loo_fits$loo_res$pred), loo_fits$loo_res_second_order$pred)
  min_scale <- min(min(loo_fits$loo_res$obs, loo_fits$loo_res$pred), loo_fits$loo_res_second_order$pred)
  
  p_stitching <- ggplot(loo_fits$loo_res, aes(x = obs, y = pred)) + geom_point() + 
    theme_bw() + geom_abline() + xlab('Observed') + ylab('Predicted') + 
    scale_x_continuous(limits = c(min_scale, max_scale)) + 
    scale_y_continuous(limits = c(min_scale, max_scale)) + 
    ggtitle('Stitching Method')
  
  p_regression <- ggplot(loo_fits$loo_res_second_order, aes(x = obs, y = pred)) + geom_point() + 
    theme_bw() + geom_abline() + xlab('Observed') + ylab('Predicted') + 
    scale_x_continuous(limits = c(min_scale, max_scale)) + 
    scale_y_continuous(limits = c(min_scale, max_scale)) + 
    ggtitle('Second-order Regression')
  
  p <- plot_grid(p_stitching, p_regression)
  
  title <- ggdraw() + 
    draw_label(
      paste0(name),
      fontface = 'bold'
      ) + 
    theme(plot.margin = margin(0, 0, 0, 7))
  p <- plot_grid(
    title, p,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  print(p)
  p_list[[dataset_id]] <- p
  
  #save(r2_loo, r2_linear_loo, res, paste0('~/global_epistasis_microbes/code/model_comparison/', name, '_model_comparison_loo.RData')) 
}


plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], p_list[[5]], p_list[[6]])

#reshape 
res_long <- res %>% pivot_longer(-name, names_to = 'r2_type') 

#set order 
res_long$r2_type <- factor(res_long$r2_type, levels = c('r2_first_order', 'r2_second_order', 'r2'))

#plot
res_long %>% 
  ggplot(aes(x = name, y = value, group = interaction(name,r2_type), fill = r2_type)) + 
  geom_col(position = 'dodge') + theme_bw() + xlab('Dataset') +
  ylab(bquote(italic(R)^2)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  guides(fill = guide_legend('Model')) +
  scale_fill_discrete(labels=c('r2_first_order' = 'First-order linear regression', 
                               'r2_second_order'= 'Second-order linear regression', 
                               'r2' = 'Stitching method'))


ggsave('~/global_epistasis_microbes/Figures/stitching_vs_reg_all_data.png', height = 8, width = 10)

################################
# run oof for all datasets #####
################################

#list all dataframe paths
path <- '~/global_epistasis_clean/Data/'
all_datasets <- list.files(path = path)

# want GE_biomass, Kuebbing natives, langenheder, pyoverdine, starch
# ghedini ~will~ throw an error. Could implement this with a more robust check
# to make sure we can fit linear patterns, but this also has issues and maybe isn't
# needed? 

all_datasets <- all_datasets[c(2,6,7,8,9)]

p_list <- list()
for (dataset_id in 1:length(all_datasets)){
  df <- read_csv(paste0(path, all_datasets[dataset_id]))
  
  name <- str_remove(all_datasets[dataset_id], '.csv')
  N <- ncol(df) - 1 
  
  res <- fit_all_models_oof(df, N, reps = 50, v = TRUE)
  
  to_plot <- res %>% dplyr::select(c(r2, r2_lin, r2_cv_pairs, infit_pct)) 
  colnames(to_plot) <- c('Stitching Method', 'First-order Linear Regression', 
                         'Second-order Linear Regression', 'infit_pct')
  to_plot <- to_plot %>% pivot_longer(-infit_pct) 
  colnames(to_plot) <- c('infit_pct', 'Model', 'r2') 
  p <- to_plot %>% ggplot(aes(x = as.factor(infit_pct), y = r2, fill = Model)) + 
    geom_boxplot() + theme_bw() + 
    xlab('Percent data in fit') + ylab(bquote(italic(R)^2)) + 
    ggtitle(name )
  
  print(p)
  
  p_list[[dataset_id]] <- p
  
}

plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], p_list[[5]] )

