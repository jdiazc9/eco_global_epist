library(tidyverse)
library(glmnet)
library(cowplot)

source('ecoFunctions.R')

##########################
####   functions    ######
###########################
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

get_all_loo_fits <- function(df, mode = 'full', v = FALSE){
  
  N <- ncol(df) -1 
  
  all_species <- colnames(df)[1:N]
  
  if (any(rowSums(df[,1:N]) == 0)){
    df <- df[-which(rowSums(df[,1:N]) == 0),]
  }
  
  communities <- sapply(1:nrow(df),
                        FUN = function(i) paste(all_species[df[i, -ncol(df)] == 1], collapse = ','))
  unique_communities <- unique(communities)
  df <- df %>% mutate(community = communities) 
  
  #add mean fitness to df 
  df <- df %>% group_by(community) %>% mutate(mean_fitness = mean(`function`)) %>% ungroup() 
  
  mean_df <- df %>% distinct(community, .keep_all = TRUE) %>% dplyr::select(-c(`function`))
  colnames(mean_df)[which(colnames(mean_df) == 'mean_fitness')] <- 'function'
  
  
  ####################
  #### stitching #####
  ####################
  
  loo_res <- tibble()
  for (exp in unique_communities){
    print(exp)
    exp_inds <- which(mean_df$community == exp)
    
    data <- mean_df[-exp_inds,] %>% dplyr::select(-community) 
    
    obsF <- mean_df[exp_inds[1],] %>% dplyr::select(-community) 
    target <- matrix2string(obsF)
    
    data <- matrix2string(data)
    ge_data <- makeGEdata(data)
    
    if (mode == 'full'){
      eps <- try(inferAllResiduals(ge_data))
      predF <- predictF_fullClosure(target$community, data, eps)
      
    } else if (mode == 'base'){
      
      predF <- predictF_base(target$community, data)
      
    }
    
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
    y <- as.matrix(df[-exp_inds,]$`function`)   
    f <- as.formula(y ~ .)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$`function`)
    f_test <- as.formula(y_test ~ .)
    x_test <- model.matrix(f_test, df[exp_inds[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
    
    tmp <- list(obs = df[exp_inds[1],]$mean_fitness, pred = y_pred_cv)
    loo_res_first_order <- rbind(loo_res_first_order, tmp)
  }
  
  #####################################################
  # first order linear regression, ridge regularization #
  #####################################################
  ridge_loo_res_first_order <- tibble()
  for (exp in unique_communities){
    exp_inds <- which(df$community == exp)
    #set up regression
    y <- as.matrix(df[-exp_inds,]$`function`)   
    f <- as.formula(y ~ .)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids, alpha = 0) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$`function`)
    f_test <- as.formula(y_test ~ .)
    x_test <- model.matrix(f_test, df[exp_inds[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
    
    tmp <- list(obs = df[exp_inds[1],]$mean_fitness, pred = y_pred_cv)
    ridge_loo_res_first_order <- rbind(ridge_loo_res_first_order, tmp)
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
    y <- as.matrix(df[-exp_inds,]$`function`)   
    f <- as.formula(y ~ .*.)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$`function`)
    f_test <- as.formula(y_test ~ .*.)
    x_test <- model.matrix(f_test, df[exp_inds[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
    
    tmp <- list(obs = df[exp_inds[1],]$mean_fitness, pred = y_pred_cv)
    loo_res_second_order <- rbind(loo_res_second_order, tmp)
  }
  
  ###############################################################
  # second order linear regression, ridge regularization #
  ###############################################################
  count <- 1
  ridge_loo_res_second_order <- tibble()
  for (exp in unique_communities){
    if (v){
      print(paste0(count , ' out of ', length(unique_communities)))
      count <- count + 1
    }
    
    exp_inds <- which(df$community == exp)
    #set up regression
    y <- as.matrix(df[-exp_inds,]$`function`)   
    f <- as.formula(y ~ .*.)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids, alpha = 0) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$`function`)
    f_test <- as.formula(y_test ~ .*.)
    x_test <- model.matrix(f_test, df[exp_inds[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
    
    tmp <- list(obs = df[exp_inds[1],]$mean_fitness, pred = y_pred_cv)
    ridge_loo_res_second_order <- rbind(ridge_loo_res_second_order, tmp)
  }
  
  loo_res$pred[loo_res$pred < 0] <- 0
  loo_res_first_order$pred[loo_res_first_order$pred < 0] <- 0
  loo_res_second_order$pred[loo_res_second_order$pred < 0] <- 0
  ridge_loo_res_first_order$pred[ridge_loo_res_first_order$pred < 0] <- 0
  ridge_loo_res_second_order$pred[ridge_loo_res_second_order$pred < 0] <- 0
  
  ### aggregate loo R^2 ###
  r2 <- cor(loo_res$obs, loo_res$pred)^2
  r2_first_order <- cor(loo_res_first_order$obs, loo_res_first_order$pred)^2
  r2_second_order <- cor(loo_res_second_order$obs, loo_res_second_order$pred)^2
  r2_ridge_first_order <- cor(ridge_loo_res_first_order$obs, ridge_loo_res_first_order$pred)^2
  r2_ridge_second_order <- cor(ridge_loo_res_second_order$obs, ridge_loo_res_second_order$pred)^2
  
  return(list(r2 = r2,   #stitching 
              r2_first_order = r2_first_order,  #regression, l1 regularization
              r2_second_order = r2_second_order, 
              r2_ridge_first_order = r2_ridge_first_order, #regression, l2 regularization
              r2_ridge_second_order = r2_ridge_second_order,
              loo_res = loo_res,  #stitching, obs vs. pred
              loo_res_first_order = loo_res_first_order, #regression, l1 regularization, obs vs.pred
              loo_res_second_order = loo_res_second_order,
              ridge_loo_res_first_order = ridge_loo_res_first_order, #regression, l2 regularization, obs vs. pred
              ridge_loo_res_second_order = ridge_loo_res_second_order))
}

# get unregularized regression R^2
get_unregularized_fits <- function(df){
  N <- ncol(df) -1 
  
  all_species <- colnames(df)[1:N]
  
  if (any(rowSums(df[,1:N]) == 0)){
    df <- df[-which(rowSums(df[,1:N]) == 0),]
  }
  
  communities <- sapply(1:nrow(df),
                        FUN = function(i) paste(all_species[df[i, -ncol(df)] == 1], collapse = ','))
  unique_communities <- unique(communities)
  df <- df %>% mutate(community = communities) 
  
  #add mean fitness to df 
  df <- df %>% group_by(community) %>% mutate(mean_fitness = mean(`function`)) %>% ungroup() 
  
  mean_df <- df %>% distinct(community, .keep_all = TRUE) %>% dplyr::select(-c(`function`))
  colnames(mean_df)[which(colnames(mean_df) == 'mean_fitness')] <- 'function'
  
  ##############################
  ## unregularized regression ##
  ##############################
  #first order
  lin_fit <- lm(`function` ~ ., df[, 1:(N + 1)])
  distinct_preds <- cbind(df, lin_preds = lin_fit$fitted.values) %>% distinct(community, .keep_all = TRUE)
  
  r2_ur_first_order <- cor(distinct_preds$mean_fitness, distinct_preds$lin_preds)^2
  
  #second order
  lin_fit_second_order <- lm(`function` ~ .*. , df[, 1:(N + 1)])
  distinct_preds <- cbind(df, lin_preds = lin_fit_second_order$fitted.values) %>% distinct(community, .keep_all = TRUE)
  
  r2_ur_second_order <- cor(distinct_preds$mean_fitness, distinct_preds$lin_preds)^2
  
  return(list(r2_ur_first_order = r2_ur_first_order, 
              r2_ur_second_order = r2_ur_second_order))
}


#################################
### get leave-one-out fits for ##
########  all datasets  #########
#################################

#list all dataframe paths
path <- '../data_sets/'
all_datasets <- list.files(path = path)
all_datasets <- c(all_datasets, all_datasets[grepl('Kuebbing', all_datasets)]) #duplicate Kuebbing
all_datasets <- c(paste0(path, all_datasets), '../pyoverdine_data/training_set.csv' )

all_datasets <- all_datasets[!grepl('Clark', all_datasets)]


#kuebbing flag 
sub_dataset1 <- 'TRUE'

p_list <- list()
res <- tibble()

lm_res <- tibble() 
res_full <- tibble()

for (dataset_id in 1:length(all_datasets)){
  df <- read_csv(all_datasets[dataset_id], show_col_types = F)
  
  name <- str_remove(basename(all_datasets[dataset_id]), '.csv')

  mode <- 'full'
  
  if (grepl('training_set', name)){
    df$`function` <- rowMeans(df[, 9:11])
    df <- df[, c(1:8, 12)]
  }
  
  if (grepl('Ghedini', name)){
    df$`function` <- df$`function`/1e4
  }
  
  if (name == 'butyrate_Clark2021'){
    mode <- 'base'
  } else {
    mode <- 'full'
  }
  
  if (name == "plant-biommass_Kuebbing2016_all"){
    if (sub_dataset1){
      df <- df[,c(1:4, 9)]
      name <- paste0(name, '_natives')
      sub_dataset1 <- FALSE
    } else{
      df <- df[,c(5:8, 9)]
      name <- paste0(name, '_invasives')
    }
  } 
  
  df <- aggregate(`function` ~ .,
                  data = df,
                  FUN = mean)
  
  loo_fits <- get_all_loo_fits(df, mode)

  
  res_full <- rbind(res_full,
                    rbind(data.frame(dataset = name,
                                     method = 'stitching',
                                     pred = loo_fits$loo_res$pred,
                                     obs = loo_fits$loo_res$obs,
                                     sq_err = abs(loo_fits$loo_res$obs - loo_fits$loo_res$pred)/mean(loo_fits$loo_res$obs)),
                          data.frame(dataset = name,
                                     method = 'first_order',
                                     pred = loo_fits$loo_res_first_order$pred,
                                     obs = loo_fits$loo_res_first_order$obs,
                                     sq_err = abs(loo_fits$loo_res_first_order$obs - loo_fits$loo_res_first_order$pred)/mean(loo_fits$loo_res_first_order$obs)),
                          data.frame(dataset = name,
                                     method = 'second_order',
                                     pred = loo_fits$loo_res_second_order$pred,
                                     obs = loo_fits$loo_res_second_order$obs,
                                     sq_err = abs(loo_fits$loo_res_second_order$obs - loo_fits$loo_res_second_order$pred)/mean(loo_fits$loo_res_second_order$obs))), 
                          data.frame(dataset = name,
                                     method = 'first_order_ridge',
                                     pred = loo_fits$ridge_loo_res_first_order$pred,
                                     obs = loo_fits$ridge_loo_res_first_order$obs,
                                     sq_err = abs(loo_fits$ridge_loo_res_first_order$obs - loo_fits$ridge_loo_res_first_order$pred)/mean(loo_fits$ridge_loo_res_first_order$obs)),
                          data.frame(dataset = name,
                                     method = 'second_order_ridge',
                                     pred = loo_fits$ridge_loo_res_second_order$pred,
                                     obs = loo_fits$ridge_loo_res_second_order$obs,
                                     sq_err = abs(loo_fits$ridge_loo_res_second_order$obs - loo_fits$ridge_loo_res_second_order$pred)/mean(loo_fits$ridge_loo_res_second_order$obs)))
  
  #save unregularized regression fits
  ur_fits <- get_unregularized_fits(df)
  lm_res <- rbind(lm_res, 
                  rbind(data.frame(dataset = name, 
                                   r2_ur_first_order = ur_fits$r2_ur_first_order, 
                                   r2_ur_second_order = ur_fits$r2_ur_second_order)))

}

legacy <- res_full

res_full$invasives <- grepl('_invasives', res_full$dataset)
res_full$dataset <- gsub('_invasives', '', res_full$dataset)
res_full$dataset <- gsub('_natives', '', res_full$dataset)

# plot pred. vs obs.
res_full$dataset <- setNames(c('Bacterial\nstarch hydrolysis',
                               'Phytoplankton\nbiomass',
                               'Above-ground\nplant biomass',
                               'Bacterial\nxylose oxidation',
                               'Bacterial\npyoverdine secretion'),
                             unique(res_full$dataset))[res_full$dataset]
res_full$dataset <- factor(res_full$dataset,
                           levels = c('Bacterial\npyoverdine secretion',
                                      'Above-ground\nplant biomass',
                                      'Phytoplankton\nbiomass',
                                      'Bacterial\nxylose oxidation',
                                      'Bacterial\nstarch hydrolysis'))

scales_limits_pred <- rbind(aggregate(pred~method+dataset,
                                      data = res_full,
                                      FUN = max),
                            aggregate(pred ~ method+dataset,
                                      data = res_full,
                                      FUN = min))
scales_limits_obs <- rbind(aggregate(obs~method+dataset,
                                     data = res_full,
                                     FUN = max),
                           aggregate(obs~method+dataset,
                                     data = res_full,
                                     FUN = min))

colnames(scales_limits_obs) <- c('method', 'dataset', 'f')
colnames(scales_limits_pred) <- c('method', 'dataset', 'f')

scales_limits <- rbind(scales_limits_pred, scales_limits_obs)

scales_limits <- rbind(aggregate(f~dataset,
                                 data = scales_limits,
                                 FUN = max),
                       aggregate(f~dataset,
                                 data = scales_limits,
                                 FUN = min))

scales_limits <- rbind(cbind(scales_limits, method = 'stitching'),
                       cbind(scales_limits, method = 'first_order'),
                       cbind(scales_limits, method = 'second_order'))

scales_limits$pred <- scales_limits$f
scales_limits$obs <- scales_limits$f

scales_limits <- scales_limits[, c('method', 'dataset', 'pred', 'obs')]
scales_limits$invasives <- FALSE

######################################################
### plot predicted versus observed values,         ###
### stitching & regression with l1 regularization  ###
######################################################

res_full %>% filter(!method %in% c('first_order_ridge', 'second_order_ridge')) %>% 
  ggplot(aes(x = obs, y = pred, shape = invasives)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_point(cex = 3) +
  geom_blank(data = scales_limits, aes(x = pred, y = obs)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = expression(paste('Observed ', italic('F'), ' [a.u.]', sep = ''))) +
  scale_y_continuous(breaks = pretty_breaks(n = 3),
                     name = expression(paste('Predicted ', italic('F'), ' [a.u.]', sep = ''))) +
  scale_shape_manual(values = c(16, 1)) +
  facet_wrap(method~dataset,
             scales = 'free',
             nrow = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

# ggsave(filename = '../plots/stitching_vs_regression.pdf',
#        device = 'pdf',
#        dpi = 600,
#        width = 300,
#        height = 200,
#        units = 'mm',
#        limitsize = F)

#############################
## plot R^2 for each model ##
#############################
get_r2 <- function(observed, predicted){
  
  RSS <- sum((observed - predicted)^2)
  TSS <- sum((observed - mean(observed))^2)
  
  return (1 - RSS/TSS)
  
}

res_full$inter <- paste(res_full$dataset, res_full$invasives, sep = ' / ')

rsq <- do.call(rbind,
               lapply(unique(res_full$inter),
                      FUN = function(inter) {
                        
                        rsq <- data.frame(dataset = character(0),
                                          method = character(0),
                                          r2 = numeric(0))
                        
                        for (method in unique(res_full$method)) {
                          #lmod <- lm(formula = pred~obs,
                          #           data = res_full[res_full$inter == inter & res_full$method == method, ])
                          res_full_tmp <- res_full[res_full$inter == inter & res_full$method == method, ]
                          lmod <- get_r2(res_full_tmp$obs, res_full_tmp$pred)
                          rsq <- rbind(rsq,
                                       data.frame(inter = inter,
                                                  method = method,
                                                  r2 = lmod)) #summary(lmod)$r.squared))
                        }
                        
                        return(rsq)
                        
                      }))

rsq$invasives <- grepl('TRUE', rsq$inter)

rsq$inter <- gsub(' / TRUE', '', rsq$inter)
rsq$inter <- gsub(' / FALSE', '', rsq$inter)
colnames(rsq)[1] <- 'dataset'

rsq$dataset <- factor(rsq$dataset,
                      levels = c('Bacterial\npyoverdine secretion',
                                 'Above-ground\nplant biomass',
                                 'Phytoplankton\nbiomass',
                                 'Bacterial\nxylose oxidation',
                                 'Bacterial\nstarch hydrolysis'))

rsq$natives <- !rsq$invasives


### R^2 values, 
### stitching method & regression with l1 regularization

rsq %>% filter(!method %in% c('first_order_ridge', 'second_order_ridge')) %>% 
  ggplot(aes(x = method, y = r2, fill = interaction(method, natives), color = interaction(method, natives))) +
  geom_bar(stat = 'identity',
           width = 0.6,
           position = 'dodge') +
  facet_wrap(~dataset,
             nrow = 1) +
  scale_y_continuous(name = 'R2',
                     limits = c(0, 1),
                     breaks = pretty_breaks(n = 4)) +
  scale_color_manual(values = c('#99d7dc', '#176766', '#b33a3b',
                                'white', 'white', 'white')) +
  scale_fill_manual(values = c('white', 'white', 'white',
                               '#99d7dc', '#176766', '#b33a3b')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1.6,
        #axis.text = element_text(size = 16),
        #axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'bottom') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.5) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5)


# ggsave(filename = '../plots/stitching_vs_regression_R2.pdf',
#        device = 'pdf',
#        dpi = 600,
#        width = 200,
#        height = 150,
#        units = 'mm',
#        limitsize = F)

##############################################
##  plot R^2 values,                        ##
##  regression with l1 & l2 regularization  ##
##############################################

rsq %>% filter(!method %in% c('stitching')) %>% 
  ggplot(aes(x = method, y = r2, fill = interaction(method, natives), color = interaction(method, natives))) +
  geom_bar(stat = 'identity',
           width = 0.6,
           position = 'dodge') +
  facet_wrap(~dataset,
             nrow = 1) +
  scale_y_continuous(name = 'R2',
                     limits = c(0, 1),
                     breaks = pretty_breaks(n = 4)) +
  scale_color_manual(values = c('#99d7dc', '#176766', '#b33a3b', 'darkred',
                                'white', 'white', 'white', 'white')) +
  scale_fill_manual(values = c('white', 'white', 'white', 'white',
                               '#99d7dc', '#176766', '#b33a3b', 'darkred')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1.6,
        #axis.text = element_text(size = 16),
        #axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'bottom') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.5) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5)

# ggsave(filename = '../plots/ridge_vs_lasso_R2.pdf',
#        device = 'pdf',
#        dpi = 600,
#        width = 200,
#        height = 150,
#        units = 'mm',
#        limitsize = F)

#########################################
##  plot unregularized regression R^2  ##
########################################
colnames(lm_res) <- c('Dataset', 'First Order', 'Second Order')

lm_res %>% pivot_longer(-Dataset) %>% ggplot(aes(x = name, y = value, fill = name)) +
  geom_bar(stat = 'identity', width = 0.6, position = 'dodge') + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1.6,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = c('gray', 'black')) +
  facet_wrap(~Dataset, nrow = 1) +
  ylab(expression(R^2)) + 
  labs(fill = '') + 
  ggtitle('Summary of unregularized regression fits')

ggsave(filename = '../plots/unregR2.pdf',
       device = 'pdf',
       dpi = 600,
       width = 230,
       height = 150,
       units = 'mm',
       limitsize = F)


# top and bottom-performing communities
threshold <- 0.1
top_bottom <- do.call(data.frame, aggregate(obs ~ inter,
                                  data = res_full,
                                  FUN = function(x) quantile(x, probs = c(threshold, 1 - threshold))))
colnames(top_bottom)[2:3] <- c('lower_bound', 'upper_bound')
which_top_bottom <- sapply(1:nrow(res_full),
                           FUN = function(i) {
                             res_full$obs[i] < top_bottom$lower_bound[top_bottom$inter == res_full$inter[i]] | res_full$obs[i] > top_bottom$upper_bound[top_bottom$inter == res_full$inter[i]]
                           })

res_top_bottom <- res_full[which_top_bottom, ]

res_top_bottom$err <- abs(res_top_bottom$pred - res_top_bottom$obs)
res_top_bottom$rel_err <- log10(abs(res_top_bottom$pred - res_top_bottom$obs)/res_top_bottom$obs)

ggplot(res_top_bottom, aes(x = method, y = rel_err, fill = interaction(method, !invasives), color = interaction(method, !invasives))) +
  geom_boxplot(outlier.colour = 'gray',
               width = 0.5) +
  facet_wrap(~dataset,
             nrow = 1,
             scales = 'free_y') +
  scale_y_continuous(name = expression(log[10]~'|'*italic(F)['pred'] - italic(F)['obs']*'|'*'/'*italic(F)['obs']),
                     breaks = pretty_breaks(n = 4)) +
  scale_color_manual(values = c('#99d7dc', '#176766', '#b33a3b',
                                'black', 'black', 'black')) +
  scale_fill_manual(values = c('white', 'white', 'white',
                               '#99d7dc', '#176766', '#b33a3b')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1.6,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'bottom') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.5) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/stitching_vs_regression_top_bottom_relErr.pdf',
       device = 'pdf',
       dpi = 600,
       width = 230,
       height = 150,
       units = 'mm',
       limitsize = F)

ggplot(res_top_bottom, aes(x = method, y = err, fill = interaction(method, !invasives), color = interaction(method, !invasives))) +
  geom_boxplot(outlier.colour = 'gray',
               width = 0.5) +
  facet_wrap(~dataset,
             nrow = 1,
             scales = 'free_y') +
  scale_y_continuous(name = expression(log[10]~'|'*italic(F)['pred'] - italic(F)['obs']*'|'*'/'*italic(F)['obs']),
                     breaks = pretty_breaks(n = 4)) +
  scale_color_manual(values = c('#99d7dc', '#176766', '#b33a3b',
                                'black', 'black', 'black')) +
  scale_fill_manual(values = c('white', 'white', 'white',
                               '#99d7dc', '#176766', '#b33a3b')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1.6,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'bottom') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.5) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/stitching_vs_regression_top_bottom_absErr.pdf',
       device = 'pdf',
       dpi = 600,
       width = 230,
       height = 150,
       units = 'mm',
       limitsize = F)

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
  
  #save(r2_loo, r2_linear_loo, res, 
    #paste0('../', name, '_model_comparison_loo.RData')) 


plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], p_list[[5]]) 

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


##ggsave('~/global_epistasis_microbes/Figures/stitching_vs_reg_all_data.png', height = 8, width = 10)
