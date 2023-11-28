library(tidyverse)
library(glmnet)
library(cowplot)

# make matrix of communities for a given number of species
makeComms <- function(N = 6, name_prefix = 'sp_') {
  
  comms <- as.data.frame(do.call(rbind,
                                 lapply(0:N,
                                        FUN = function(i) {
                                          
                                          sp <- t(combn(N, i))
                                          
                                          comms_i <- do.call(rbind,
                                                             lapply(1:nrow(sp),
                                                                    FUN = function(i) {
                                                                      
                                                                      comm_i <- rep(0, N)
                                                                      comm_i[sp[i, ]] <- 1
                                                                      return(comm_i)
                                                                      
                                                                    }))
                                          
                                          return(comms_i)
                                          
                                        })))
  colnames(comms) <- paste(name_prefix, 1:N, sep = '')
  
  return(comms)
  
}

# make landscape from given coefficients of Fourier (or Taylor) representation
funFromCoeff <- function(comms, coeff, base = 'Fourier') {
  
  fun <- sapply(1:nrow(comms),
                FUN = function(i) {
                  
                  community <- comms[i, ]
                  if (base == 'Fourier') community[community == 0] <- -1 # Fourier representation (x = +1/-1), default input is in Taylor representation (x = 0/1)
                  sigma <- unlist(lapply(0:length(community),
                                         FUN = function(i) apply(t(combn(ncol(comms), i)),
                                                                 FUN = function(x) prod(community[x]),
                                                                 MARGIN = 1)))
                  return(sum(sigma*coeff))
                  
                })
  
  return(fun)
  
}

# fraction of functional variance explained by coefficients of each order
varFractionFromCoeff <- function(coeff) {
  
  total_var <- sum(coeff[2:length(coeff)]^2)
  var_by_order <- aggregate(value ~ order,
                            data = data.frame(order = nchar(names(coeff[2:length(coeff)])), value = coeff[2:length(coeff)]),
                            FUN = function(x) sum(x^2))$value / total_var
  
  return(setNames(var_by_order, 1:length(var_by_order)))
  
}


# get distribution of epistatic (and fitness) effects
getDEE <- function(df) {
  
  N <- ncol(df) - 1
  species <- colnames(df)[1:N]
  pairs <- t(combn(N, 2))
  pairs <- rbind(pairs, cbind(pairs[, 2], pairs[, 1]))
  
  DEE <- do.call(rbind,
                 lapply(1:nrow(pairs),
                        FUN = function(p) {
                          
                          out <- data.frame(species_i = species[pairs[p, 1]],
                                            species_j = species[pairs[p, 2]],
                                            background = NA,
                                            fun_background = NA,
                                            dF_i = NA,
                                            dF_j = NA,
                                            eps_ij = NA,
                                            fun_ij = NA)
                          
                          C_0 <- df[df[, pairs[p, 1]] == 0 & df[, pairs[p, 2]] == 0, ]
                          C_i <- df[df[, pairs[p, 1]] == 1 & df[, pairs[p, 2]] == 0, ]
                          C_j <- df[df[, pairs[p, 1]] == 0 & df[, pairs[p, 2]] == 1, ]
                          C_ij <- df[df[, pairs[p, 1]] == 1 & df[, pairs[p, 2]] == 1, ]
                          
                          if (min(c(nrow(C_0), nrow(C_i), nrow(C_j), nrow(C_ij)))) {
                            
                            colnames(C_ij)[ncol(C_ij)] <- 'fun_ij'
                            C_ij[, pairs[p, 1]] <- 0
                            C_ij[, pairs[p, 2]] <- 0
                            
                            
                            colnames(C_i)[ncol(C_i)] <- 'fun_i'
                            C_i[, pairs[p, 1]] <- 0
                            
                            
                            colnames(C_j)[ncol(C_j)] <- 'fun_j'
                            C_j[, pairs[p, 2]] <- 0
                            
                            
                            colnames(C_0)[ncol(C_0)] <- 'fun_0'
                            
                            C_all <- merge(merge(merge(C_0, C_i), C_j), C_ij)
                            dF_i <- C_all$fun_i - C_all$fun_0
                            dF_j <- C_all$fun_j - C_all$fun_0
                            eps_ij <- C_all$fun_ij - C_all$fun_i - C_all$fun_j + C_all$fun_0

                            if (nrow(C_all)) out <- data.frame(species_i = species[pairs[p, 1]],
                                                               species_j = species[pairs[p, 2]],
                                                               background = matrix2string(C_all[, 1:(N+1)])$community,
                                                               fun_background = C_all$fun_0,
                                                               dF_i = dF_i,
                                                               dF_j = dF_j,
                                                               eps_ij = eps_ij,
                                                               fun_ij = C_all$fun_ij)
                            
                          }
                          
                          return(out)
                          
                        }))
  
  return(DEE)
  
}

# wrapper function: evaluate performance of predictive method via leave-one-out cross-validation
evaluatePredictions <- function(df) {
  
  pred_vs_obs <- do.call(rbind,
                         lapply(2:nrow(df),
                                FUN = function(n_oos) {
                                  
                                  df_oos <- df[n_oos, , drop = F]
                                  df_is <- df[-n_oos, , drop = F]
                                  
                                  gedf <- makeGEdata(matrix2string(df_is))
                                  theta <- inferAllResiduals(gedf)
                                  
                                  out <- merge(matrix2string(df_oos),
                                               predictF_fullClosure(matrix2string(df_oos)$community, matrix2string(df_is), theta),
                                               by = 'community',
                                               suffixes = c('_true', '_predicted'))
                                  return(out)
                                  
                                }))
  
  return(pred_vs_obs)
  
}

# same as before but without residual inference
evaluatePredictions_base <- function(df) {
  
  pred_vs_obs <- do.call(rbind,
                         lapply(2:nrow(df),
                                FUN = function(n_oos) {
                                  
                                  df_oos <- df[n_oos, , drop = F]
                                  df_is <- df[-n_oos, , drop = F]
                                  
                                  gedf <- makeGEdata(matrix2string(df_is))
                                  
                                  out <- merge(matrix2string(df_oos),
                                               predictF_base(matrix2string(df_oos)$community, matrix2string(df_is)),
                                               by = 'community',
                                               suffixes = c('_true', '_predicted'))
                                  return(out)
                                  
                                }))
  
  return(pred_vs_obs)
  
}

# wrapper function: evaluate performance of predictive method leaving out a fraction of the sample
evaluatePredictions_mult <- function(df, fraction_out_of_sample = 0.6) {
  
  which_empty <- sapply(1:nrow(df),
                        FUN = function(i) all(df[i, 1:(ncol(df) - 1)] == 0))
  comm_index <- 1:nrow(df)
  comm_index <- comm_index[!which_empty]
  comm_insample <- sort(c(which(which_empty),
                        sample(comm_index, size = floor(nrow(df)*(1 - fraction_out_of_sample)))))
  
  df_is <- df[comm_insample, ]
  df_oos <- df[-comm_insample, ]
  
  gedf <- makeGEdata(matrix2string(df_is))
  theta <- inferAllResiduals(gedf)
  
  out <- merge(matrix2string(df_oos),
               predictF_fullClosure(matrix2string(df_oos)$community, matrix2string(df_is), theta),
               by = 'community',
               suffixes = c('_true', '_predicted'))
  
  return(out)
  
}

# wrapper function: plot FEEs, clean
plotFEEs_clean <- function(df) {
  
  gedf <- makeGEdata(matrix2string(df))
  fees <- makeFEEs(gedf)
  fees <- cbind(knock_in = rownames(fees),
                fees)
  
  myplot <- 
    ggplot(gedf,
           aes(x = background_f, y = d_f)) +
    geom_abline(slope = 0, intercept = 0, color = 'gray') +
    geom_point() +
    geom_abline(data = fees,
                aes(slope = b, intercept = a, color = b)) +
    scale_color_gradient2(low = 'firebrick1',
                          high = 'deepskyblue',
                          mid = 'black') +
    scale_x_continuous(name = 'Background function [a.u.]') +
    scale_y_continuous(name = expression(paste(Delta*italic(F), ' [a.u.]'))) +
    facet_wrap(~ knock_in,
               ncol = 2,
               labeller = labeller(knock_in = setNames(paste('species ', 1:N, sep = ''),
                                                       paste('sp_', 1:N, sep = '')))) +
    theme_bw() +
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 14),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  print(myplot)
  
}

# wrapper function to plot a fitness graph
plotFitnessGraph <- function(df) {

  # make edges of fitness graph
  genots <- matrix2string(df)$community
  
  # number of mutations in genotype
  nMut <- function(genot) sapply(genot,
                                 FUN = function(genot_i) length(strsplit(genot_i, split = ',')[[1]]))
  
  # check if a genotype is descendant from another
  isDescendant <- function(this_genot, of_this_genot) {
    
    this_genot <- strsplit(this_genot, split = ',')[[1]]
    of_this_genot <- strsplit(of_this_genot, split = ',')[[1]]
    
    return(all(of_this_genot %in% this_genot))
    
  }
  
  #make edges
  makeEdges <- function(genots) {
    
    edges <- data.frame(source = character(0),
                        target = character(0),
                        source.nmut = numeric(0),
                        target.nmut = numeric(0))
    
    for(s in genots) {
      
      t <- genots[sapply(genots,
                         isDescendant,
                         of_this_genot = s) & nMut(genots) == nMut(s)+1]
      if(length(t)) {
        edges <- rbind(edges,
                       data.frame(source = s,
                                  target = as.character(t),
                                  source.nmut = as.numeric(nMut(s)),
                                  target.nmut = as.numeric(nMut(s)) + 1))
      }
      
    }
    
    edges <- cbind(edge_id = paste('edge_', 1:nrow(edges), sep = ''),
                   edges)
    
    return(edges)
  
  }
  
  # plot landscape
  plotGraph <- function(df, save.plot = F) {
    
    mycolors <- c('#939598', '#d68f28', '#415ba9', '#a96cad')
    
    n_mut <- ncol(df) - 1
    
    landscape <- matrix2string(df)
    colnames(landscape) <- c('genot', 'f')
    
    edges <- makeEdges(landscape$genot)
    
    df <- cbind(edges,
                source.f = setNames(landscape$f, landscape$genot)[edges$source],
                target.f = setNames(landscape$f, landscape$genot)[edges$target])
    df$source.f[is.na(df$source.f)] <- landscape$f[landscape$genot == '']
    
    if ('color' %in% colnames(landscape)) {
      df <- merge(df, landscape[, c('genot', 'color')], by.x = 'target', by.y = 'genot')
    } else {
      df$color <- 'A'
    }
    df <- df[, c('edge_id', 'source', 'target', 'source.nmut', 'target.nmut', 'source.f', 'target.f', 'color')]
    
    dfx <- gather(df[, c(1, 4, 5)], position, nmut, source.nmut:target.nmut)
    dfx$position <- setNames(c('source', 'target'), c('source.nmut', 'target.nmut'))[dfx$position]
    
    dfy <- gather(df[, c(1, 6, 7)], position, f, source.f:target.f)
    dfy$position <- setNames(c('source', 'target'), c('source.f', 'target.f'))[dfy$position]
    
    dfxy <- merge(dfx, dfy, by = c('edge_id', 'position'))
    
    df <- merge(dfxy, df[, c('edge_id', 'color')], by = 'edge_id')
    
    dy <- min(c(max(landscape$f) - landscape$f[1], landscape$f[1] - min(landscape$f)))
    dy <- round(dy/0.1)*0.1
    ybreaks <- seq(landscape$f[1] - 10*dy, landscape$f[1] + 10*dy, by = dy)
    
    myplot <-
      ggplot(df, aes(x = nmut, y = f, group = edge_id, color = color)) +
      geom_abline(slope = 0,
                  intercept = landscape$f[landscape$genot == ''],
                  color = '#d1d3d4') +
      geom_line() +
      scale_x_continuous(name = '# of species',
                         breaks = 0:n_mut,
                         labels = as.character(0:n_mut)) +
      scale_y_continuous(name = 'Function [a.u.]',
                         breaks = pretty_breaks(n = 3),
                         expand = c(0.05, 0.05)) +
      scale_color_manual(values = setNames(mycolors, LETTERS[1:length(mycolors)])) +
      theme_bw() +
      theme(aspect.ratio = 0.6,
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = 'none',
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16)) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
    
    if (save.plot != F) {
      ggsave(myplot,
             file = paste('../plots/', save.plot, '.pdf', sep = ''),
             dpi = 600,
             width = 100,
             height = 80,
             units = 'mm')
    }
    
    return(myplot)
    
  }
  
  return(plotGraph(df))

}

# inverse logistic function
myLogistic <- function(x, A, B) {
  
  y <- (0.5 - 0.5*tanh(A*(x - B)))
  
  #y <- 1 / (1 + exp(A*(x - B)^1))
  
  return(y)
}


### ABBY'S WRAPPER FUNCTIONS

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
  
  # this next line removes the 0 community 
  # so as not to predict it out of fit & 
  # keep an appropriate comparison to the stitching method
  
  if (any(rowSums(df[,1:N]) == 0)){
    df <- df[-which(rowSums(df[,1:N]) == 0),]
  }
  
  communities <- sapply(1:nrow(df),
                        FUN = function(i) paste(all_species[df[i, -ncol(df)] == 1], collapse = ','))
  unique_communities <- unique(communities)
  df <- df %>% mutate(community = communities) 
  
  #add mean fitness to df 
  df <- df %>% group_by(community) %>% mutate(mean_fitness = mean(`fun`)) %>% ungroup() 
  
  mean_df <- df %>% distinct(community, .keep_all = TRUE) %>% dplyr::select(-c(`fun`))
  colnames(mean_df)[which(colnames(mean_df) == 'mean_fitness')] <- 'function'
  
  #################################
  # first order linear regression #
  #################################
  loo_res_first_order <- tibble()
  for (exp in unique_communities){
    exp_inds <- which(df$community == exp)
    
    #set up regression
    y <- as.matrix(df[-exp_inds,]$`fun`)   
    f <- as.formula(y ~ .)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$`fun`)
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
    y <- as.matrix(df[-exp_inds,]$`fun`)   
    f <- as.formula(y ~ .*.)
    x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
    
    unique_communities_train <- unique_communities[-which(unique_communities == exp)]
    fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
    
    cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
    
    #get out of fit cv 
    y_test <- as.matrix(df[exp_inds[1],]$`fun`)
    f_test <- as.formula(y_test ~ .*.)
    x_test <- model.matrix(f_test, df[exp_inds[1],] %>% dplyr::select(all_of(all_species)))
    y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
    
    tmp <- list(obs = df[exp_inds[1],]$mean_fitness, pred = y_pred_cv)
    loo_res_second_order <- rbind(loo_res_second_order, tmp)
  }
  
  r2_first_order <- cor(loo_res_first_order$obs, loo_res_first_order$pred)^2
  r2_second_order <- cor(loo_res_second_order$obs, loo_res_second_order$pred)^2
  
  return(data.frame(community = unique_communities,
                    fun_true = loo_res_first_order$obs,
                    fun_pred_1st = loo_res_first_order$pred,
                    fun_pred_2nd = loo_res_second_order$pred))
  
}


#######################
# calculate r/s metric #
########################

rmse <- function(observed, predicted){
  return(sqrt(mean((observed - predicted)^2)))
}


get_rs <- function(df){
  
  ## here assuming function column is named 'fitness', can be changed ##
  lin_fit <- lm(`fun` ~ . , df)
  add_coefs <- coef(lin_fit)[2:ncol(df)]
  
  res_linear <- data.frame(obs = df$`fun`, 
                           predicted = lin_fit$fitted.values)
  
  r <- rmse(res_linear$obs, res_linear$predicted)
  s <- mean(abs(add_coefs))
  
  rs <- r/s 
  
  return(rs)
}

# calculate fraction of functional variance explained by HOIs
get_vH <- function(df){
  
  lin_fit <- lm(`fun` ~ .*. , df)
  add_coefs <- coef(lin_fit)[2:ncol(df)]
  
  res_linear <- data.frame(obs = df$`fun`, 
                           predicted = lin_fit$fitted.values)
  
  vH <- 1 - var(res_linear$predicted) / var(res_linear$obs)
  
  return(vH)
}




# adapted from Abby's code: evaluate predictions of 1st and 2nd order regression leaving a fraction of observations out of sample
get_all_loo_fits_mult <- function(df, fraction_out_of_sample = 0.6, v = FALSE){
  
  N <- ncol(df) -1
  
  all_species <- colnames(df)[1:N]
  
  # this next line removes the 0 community 
  # so as not to predict it out of fit & 
  # keep an appropriate comparison to the stitching method
  
  if (any(rowSums(df[,1:N]) == 0)){
    df <- df[-which(rowSums(df[,1:N]) == 0),]
  }
  
  communities <- sapply(1:nrow(df),
                        FUN = function(i) paste(all_species[df[i, -ncol(df)] == 1], collapse = ','))
  unique_communities <- unique(communities)
  df <- df %>% mutate(community = communities) 
  
  #add mean fitness to df 
  df <- df %>% group_by(community) %>% mutate(mean_fitness = mean(`fun`)) %>% ungroup() 
  
  mean_df <- df %>% distinct(community, .keep_all = TRUE) %>% dplyr::select(-c(`fun`))
  colnames(mean_df)[which(colnames(mean_df) == 'mean_fitness')] <- 'function'
  
  #################################
  # first order linear regression #
  #################################
  loo_res_first_order <- tibble()
  
  exp <- sample(unique_communities, size = floor(fraction_out_of_sample*length(unique_communities)))
  
  exp_inds <- sort(which(sapply(df$community,
                                FUN = function(x) x %in% exp)))
    
  #set up regression
  y <- as.matrix(df[-exp_inds,]$`fun`)   
  f <- as.formula(y ~ .)
  x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
  
  unique_communities_train <- df$community[sapply(df$community,
                                                  FUN = function(x) !(x %in% names(exp_inds)))]
  fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
  
  cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
  
  #get out of fit cv 
  y_test <- as.matrix(df[exp_inds,]$`fun`)
  f_test <- as.formula(y_test ~ .)
  x_test <- model.matrix(f_test, df[exp_inds,] %>% dplyr::select(all_of(all_species)))
  y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
  
  tmp <- list(obs = df[exp_inds,]$mean_fitness, pred = y_pred_cv)
  loo_res_first_order <- rbind(loo_res_first_order, tmp)
  
  #################################
  # second order linear regression #
  #################################
  loo_res_second_order <- tibble()
    
  #set up regression
  y <- as.matrix(df[-exp_inds,]$`fun`)   
  f <- as.formula(y ~ .*.)
  x <- model.matrix(f, df[-exp_inds,] %>% dplyr::select(all_of(all_species)))
  
  fold_ids <- get_folds(df[-exp_inds,], unique_communities_train, n_folds = 10)
  
  cv_fit <- cv.glmnet(x = x, y = y, foldid = fold_ids) 
  
  #get out of fit cv 
  y_test <- as.matrix(df[exp_inds,]$`fun`)
  f_test <- as.formula(y_test ~ .*.)
  x_test <- model.matrix(f_test, df[exp_inds,] %>% dplyr::select(all_of(all_species)))
  y_pred_cv <- predict(cv_fit, x_test, s = "lambda.1se" )
  
  tmp <- list(obs = df[exp_inds,]$mean_fitness, pred = y_pred_cv)
  loo_res_second_order <- rbind(loo_res_second_order, tmp)
  
  r2_first_order <- cor(loo_res_first_order$obs, loo_res_first_order$pred)^2
  r2_second_order <- cor(loo_res_second_order$obs, loo_res_second_order$pred)^2
  
  return(data.frame(community = exp,
                    fun_true = loo_res_first_order$obs,
                    fun_pred_1st = loo_res_first_order$pred,
                    fun_pred_2nd = loo_res_second_order$pred))
  
}


