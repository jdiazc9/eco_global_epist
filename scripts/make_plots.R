#################################
# make plots 
# - coupling 
####################


mydir <- '5'
path <- "~/whynot_composition/scripts/squares/data_files/CRM_ex_supply/squares_5/"
res_path <- paste0(path, 'square_res/')
res_files <- list.files(res_path)

# if 2
#res_files <- res_files[-7]

env_ti_res_all <- tibble() 

pc_env_list <- list()
p_comp_list <- list()
p_comp_2_list <- list() 

across_mat_res <- tibble()
extinct_res <- c()

#for all files 
for (file_i in 1){ #:length(res_files)){
  print(file_i)
  
  #read in 
  load(paste0(res_path, res_files[file_i]), v  = FALSE)
  
  #check some dims 
  n <- nrow(C)
  L <- k <- ncol(C)
  
  #jacobian? 
  K <- rep(5, L)
  Cinvd <- MASS::ginv(C) %*% params$m
  xstar <- MASS::ginv(t(C)) %*% ((K - Cinvd )/Cinvd)
  J <- matrix(0, nrow = n + k, ncol = n + k)
  J[1:n, (n + 1): (n + k)] <- diag(c(xstar)) %*% C
  J[(n + 1): (n + k), (n + 1): (n + k)] <- diag(c(-K/Cinvd))
  J[(n + 1): (n + k), 1:n] <- -diag(c(Cinvd)) %*% t(C)
  
  max_jval <- max(eigen(J)$values)
  
  #might have to do this for now
  p <- plot_heatmap(C, kill_legend = TRUE) 
  
  #CTinv <- MASS::ginv(t(C)) %*% diag(1)
  #plot_heatmap(diag(1/rowSums(CTinv)) %*% CTinv)
  
  cct_norm <- C %*% t(C)  #diag(c(1/diag(C %*% t(C)))) %*% C %*% t(C) 
  plot_heatmap(C[hclust(dist(cct_norm))$order,])
  cct_norm <- cct_norm[hclust(dist(cct_norm))$order, hclust(dist(cct_norm))$order ]
  plot_heatmap(cct_norm)
  
  sp_res <- sp_out
  
  extinct <- sapply(1:nrow(sp_out), function(row) any(sp_out[row,] < THRESH))
  extinct_res[file_i] <- length(which(extinct))
  
  if (length(which(extinct) > 0)){
    sp_res <- sp_res[-which(extinct), ]
  }
  
  #plot pca spectra
  #pc_sp <- summary(prcomp(sp_res, scale = TRUE, center = TRUE))
  pc_sp <- summary(prcomp(sp_delta, scale = TRUE, center = FALSE))
  
  prot <- plot_heatmap(pc_sp$rotation, kill_legend = TRUE, my_fontsize = 6, my_title = 'PCA modes')
  
  pspec <- data.frame(pc_sp$importance[2,], s = seq(1:length(pc_sp$importance[2,]))) %>% 
    pivot_longer(-s) %>% 
    ggplot(aes(x = s, y = value)) + geom_point(size = 1) + geom_path() + theme_classic() + 
    theme(text = element_text(size = 5)) +
    xlab('Principal Component') + ylab('Variance Explained') + 
    ylim(c(0, NA)) + 
    ggtitle('Variance explained by PCs') 
  
  #pca modes heatmap vs pca spectra
  p_comp <- plot_grid(prot, pspec)
  
  ##### coupling ###
  bb <- c()
  for (my_row in 1:nrow(t(C))){
    
    bb <- c(bb, norm(t(C)[my_row,], "2"))
  }
  
  #traits
  cc <- diag(c(1/bb)) %*% t(C) %*% pc_sp$rotation
  plot(pc_sp$sdev^2, colMeans(abs(cc)))
  
  #biomass
  plot(abs(colSums(pc_sp$rotation)))
  plot(pc_sp$sdev^2, colMeans(abs(pc_sp$rotation)))
  
  #sp abundnce
  indes <- rep(1:ncol(pc_sp$rotation), ncol(pc_sp$rotation))
  plot(indes, abs(pc_sp$rotation))
  
  plot(pc_sp$sdev^2, t(abs(pc_sp$rotation))[,1])
  
  df <- data.frame(lambda = rep(pc_sp$sdev^2, 3), 
                   couplings = c(colMeans(abs(cc)), abs(colSums(pc_sp$rotation)), abs(pc_sp$rotation[2,])), #colMeans(abs(pc_sp$rotation))),
                   fn = c(rep('Traits', 10), rep('Biomass', 10), rep('Species abd', 10))) 
  df <- df %>% group_by(fn) %>% mutate(norm_couplings = couplings/sum(couplings))
  df %>% group_by(fn)%>% summarize(var(couplings))
  
  p_cc <- df %>% #filter(fn != 'Species abd') %>% 
    ggplot(aes(x = lambda, y = couplings, color = fn)) + geom_point() + 
    theme_bw() + 
    facet_wrap(~fn, scales = 'fixed') + 
    theme(legend.position = 'none') + 
    xlab(expression(lambda)) + 
    ylab('Coupling') + 
    theme(text=element_text(size= 15), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(p_cc)
  
  ggsave(p_cc, filename = '~/whynot_composition/figures/coupling_comp.pdf', width = 6, height = 4)
  
  # coupling sums to save and test 
  coups_sum <- df %>% group_by(fn) %>% mutate(coups = sum(lambda*couplings)) %>% distinct(coups)
  
  
  #plot(pc_sp$importance[2,], sort(rowMeans(sp_to_ti_res[,2:ncol(sp_to_ti_res)]), decreasing = TRUE))
  #plot((1/pc_sp$importance[1,]^2)/sum(1/pc_sp$importance[1,]^2))
  #plot(sort(rowMeans(sp_to_ti_res[,2:ncol(sp_to_ti_res)]), decreasing = FALSE))
  
  #make some intermediate traits
  inter_Ti_res <- Ti_res %*% (matrix(runif(L*L), nrow = L, ncol = L) * matrix(sample(c(0,1), L*L, replace = TRUE), nrow = L, ncol = L)) 
  env_inter_ti_res <- generate_pred_res(pert_out, inter_Ti_res, names = c('id', paste0('B', seq(1:ncol(Ti_res)))))
  sp_inter_ti_res <- generate_pred_res(sp_out, inter_Ti_res, names = c('id', paste0('B', seq(1:ncol(Ti_res)))))
  
  #biomass, here with total abundance normalization
  biomass_res <- matrix(rowSums(sp_out) )
  #biomass, here with total resource uptake normalization
  biomass_res <- matrix(rowSums(Ti_res) )
  env_biomass_res <- generate_pred_res(pert_out, biomass_res, names = c('id', paste0('bio', 1)))
  p_comp <- plot_pred_res_comp(list(env_ti_res, env_sp_res,  env_biomass_res),
                               names_list = c( 'Traits','Species',  'Biomass' ), 
                               my_fontsize = 14) 
  
  p_comp
  
  p_coup_plus_auc <- plot_grid(p_comp, p_cc)
  ggsave(p_coup_plus_auc, filename = '~/whynot_composition/figures/coup_and_auc.pdf', width = 8, height = 4)
  
  sp_to_bio_res <- generate_pred_res(sp_out, biomass_res, names = c('id', paste0('T', seq(1:ncol(biomass_res)))) )
  
  p_comp_biomass <- plot_pred_res_comp(list(env_biomass_res,  sp_to_bio_res),
                                       names_list = c( 'Environment to biomass',  'Species to biomass' ), 
                                       my_xlab = 'Number of variables used')
  
  p_comp_biomass
  
  #save(p_comp_biomass, file = paste0(save_path, 'p_comp_biomass.RData'))
  
  #plot environmental pred
  env_ti_res <- generate_pred_res(pert_out, Ti_res, names = c('id', paste0('T', seq(1:ncol(Ti_res)))))
  env_sp_res <- generate_pred_res(pert_out, sp_out, names = c('id', paste0('sp', seq(1:ncol(sp_res)))))
  
  sp_to_ti_res <- generate_pred_res(sp_out, Ti_res, names = c('id', paste0('T', seq(1:ncol(Ti_res)))) )
  
  p_comp <- plot_pred_res_comp(list(env_ti_res, env_sp_res, env_inter_ti_res),
                               names_list = c( 'Traits','Species', 'Combination traits' ))
  
  print(p_comp)
  
  p_comp_list[[file_i]] <- p_comp
  
  p_comp <- plot_pred_res_comp(list(env_ti_res, env_sp_res))
  print(p_comp)
  
  p_comp_2 <- plot_pred_res_comp(list(env_ti_res, sp_to_ti_res), names_list = c('Environment to traits', 'Species to traits'))
  print(p_comp_2)
  
  p_comp_2_list[[file_i]] <- p_comp_2
  
  ##########################################################
  # prediction from environment vs. prediction from species #
  ###########################################################
  
  p_comp_2 <- plot_pred_res_comp(list(env_ti_res, sp_to_ti_res), names_list = c('Environment to traits', 'Species to traits'))
  print(p_comp_2)
  
  p_comp_biomass <- plot_pred_res_comp(list(env_biomass_res,  sp_to_bio_res),
                                       names_list = c( 'Environment to biomass',  'Species to biomass' ), 
                                       my_xlab = 'Number of variables used')
  
  
  #### pca spectra and env res 
  p_pc_env <- plot_grid(pspec, p_comp)
  
  print(p_pc_env)
  
  #save this comparison plot for all of them
  pc_env_list[[file_i]] <- p_pc_env
  
  Cinvd <-  MASS::ginv(C) %*% params$m 
  mat <- MASS::ginv(t(C)) %*% diag(c(1/Cinvd)) #%*% diag(c(k_noise/Cinvd)^2) %*% MASS::ginv(C)
  pert_eigvals <- svd(mat)$d #eigen(mat)$values 
  pert_eigvec <- svd(mat)$u[,1] #eigen(mat)$vectors[,1]
  pert_IPR <- sum(pert_eigvec^4)
  pc_IPR <- sum(pc_sp$rotation[,1]^4)
  #pert_eigvals <- eigen(MASS::ginv(t(C)) %*% MASS::ginv(C))$values
  #pert_eigvals <- svd(MASS::ginv(t(C)) %*% diag(c(1/Cinvd^2)))$d^2
  
  cor(pert_eigvals/sum(pert_eigvals), pc_sp$importance[2,])
  
  #plot predicted vs. expected pca spectra. only works when scale, center = FALSE
  #plot((pc_sp$importance[2,]), (pert_eigvals/sum(pert_eigvals)))
  #abline(0,1)
  
  plot_heatmap(cov(sp_res))
  pert_cov_mat <- MASS::ginv(t(C)) %*% diag(c(k_noise^2/Cinvd^2)) %*% MASS::ginv(C)
  p_an <- plot_heatmap(pert_cov_mat)
  
  eigs <- pc_out$sdev^2
  #eigs[-1] <- 0
  pca1.cov <- pc_out$rotation %*% diag(eigs) %*% t(pc_out$rotation)
  p_pc <- plot_heatmap(pca1.cov)
  
  plot_grid(p_pc, p_an)
  
  plot_heatmap(pca1.cov/pert_cov_mat)
  
  mean_preds <- rowMeans(env_sp_res[, 2:(n + 1)])
  AUC <- trapz(env_sp_res$id, mean_preds)/9
  
  AUC_bio <- trapz(env_biomass_res$id, (env_biomass_res[, 2]) )/9
  AUC_sp_traits <- trapz(sp_to_ti_res$id, rowMeans(sp_to_ti_res[,2:(n + 1)]))
  
  AUC_pca <- trapz(1:n, pc_sp$importance[3,])/9
  
  this_mat <- MASS::ginv(t(C)) %*% diag(c(1/(MASS::ginv(C) %*% params$m)))
  norm_mat <- diag(c(1/rowSums(this_mat^2))) %*% this_mat^2
  mean_norm <- mean(apply(norm_mat, 1, "max"))
  
  colSums(this_mat)
  
}



#################################
# plot relative fluctuations    #
# as a function of system noise #
##################################

#### do this variance test for one system 

#load the file
file_i <- 1
load(paste0(res_path, res_files[file_i]), v  = FALSE)

params <- reset_params(C)
full_eq <- MASS::ginv(t(C)) %*% ((params$K - MASS::ginv(C) %*% params$m)/(MASS::ginv(C) %*% params$m))
survived <- which(full_eq[1:n] > .5) 

k_n_list <- c(.01, .03, .05, .07, .09)
sp_k_n_list <- tibble()
sp_plot_list <- tibble()

for (k_n in k_n_list){
  # perturb environment 
  sp_out_res <- perturb_k_square(params, start_scheme = 'uniform', 
                                 init = rep(1, params$n), 
                                 nreps = 50,
                                 k_noise = k_n, 
                                 v = TRUE, rand = TRUE)
  sp_out <- sp_out_res$sp_res
  res_out <- sp_out_res$r_res
  pert_out <- t(sp_out_res$pert_res)
  Ti_res <- sp_out %*% C
  
  all_sp_var <- c()
  for (i in 1:ncol(sp_out)){
    all_sp_var[i] <- var(sp_out[,i])/colMeans(sp_out)[i]
  }
  
  all_ti_var <- c()
  for (i in 1:ncol(Ti_res)){
    all_ti_var[i] <- var(Ti_res[,i])/colMeans(Ti_res)[i]
  }
  
  tmp <- list(type = 'species', k_noise = k_n, mean_sp_var = mean(all_sp_var), sd_sp_var = sd(all_sp_var))
  sp_plot_list <- rbind(sp_plot_list, tmp)
  
  tmp <- list(type = 'trait', k_noise = k_n, mean_sp_var = mean(all_ti_var), sd_sp_var = sd(all_ti_var))
  sp_plot_list <- rbind(sp_plot_list, tmp)
  
  #which of these guys are weak 
  colSums((sp_out < THRESH)*1) 
  extinct <- sapply(1:nrow(sp_out), function(row) any(sp_out[row,] < THRESH))
  print(paste0('extinctions', length(which(extinct))))
  
  sp_k_n_list <- rbind(sp_k_n_list, sp_out )
}


sp_plot_list %>% ggplot(aes(x = as.factor(k_noise), y = mean_sp_var, color = type)) + 
  geom_point(position=position_dodge(width = .5)) + 
  geom_linerange(aes(ymin = mean_sp_var - sd_sp_var, ymax = mean_sp_var + sd_sp_var), position=position_dodge(width = .5)) +
  xlab('Magnitude of environmental noise') + 
  ylab('Relative fluctuations') + 
  theme_classic() + 
  theme(legend.title=element_blank(), text = element_text(size = 14)) + 
  scale_color_discrete(labels = c('Species', 'Traits'))

ggsave(file = '~/whynot_composition/figures/flucs_vs_mag_noise.pdf', 
       height = 4, width = 6)

####################################
# Normalized coupling plot fig 2c  #
####################################
make_norm_coupling_plot <- function(C, pc_out, r2_df, 
                                    coup_type = 'Ti',
                                    r2_flag = TRUE, make_norm = TRUE){
  
  desired_Ti_var <- c()
  for (i in 1:ncol(Ti_res)){
    desired_Ti_var <- c(desired_Ti_var, summary(lm(Ti_res[,2] ~ sp_out %*% pc_out$rotation[,i] - 1))$r.squared)
  }
  
  pc_out <- summary(prcomp(sp_out, scale = FALSE, center = FALSE))
  summary(lm(Ti_res[,1] ~ sp_out %*% pc_out$rotation[,1]))
  
  if (coup_type == 'Ti'){
    coup_full <- abs(t(C) %*% pc_out$rotation) 
    
    coups <-  coup_full %*% diag(pc_out$importance[1,])
    coups <- apply(coups, 1, function(row) (row^2)/sum(row^2)) 
    
    coup_df <- data.frame(pc_var = pc_out$importance[2,], coups)  #t(coups))
    colnames(coup_df) <- c('pc_var', paste0('Ti_', seq(1:ncol(C))))
    
    coup_df <- coup_df %>% pivot_longer(-pc_var) 
    var_res <- apply(Ti_res, 2, function(col) var(col))
    
  } else if (coup_type == 'biomass'){
    coup_full <- abs(colSums((t(C) %*% pc_out$rotation)))
    
    coups <-  coup_full %*% diag(pc_out$importance[1,])
    coups <- apply(coups, 1, function(row) (row^2)/sum(row^2)) 
    
    coup_df <- data.frame(pc_var = pc_out$importance[2,], coups)
    colnames(coup_df) <- c('pc_var', 'biomass')
    
    coup_df <- coup_df %>% pivot_longer(-pc_var) 
  }
  
  if (r2_flag){
    r2_df <- data.frame(name = paste0('Ti_', seq(1:ncol(C))), color = r2_all_Ti_res, var_res = var_res)
    coup_df <- left_join(coup_df, r2_df, by = 'name')
  }
  
  if (r2_flag & make_norm){
    p <- coup_df %>% 
      ggplot(aes(x = pc_var, y = value, group = name, color = color)) + 
      geom_point() + 
      #geom_path() + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_color_continuous(type = 'viridis') + 
      xlab('Percent compositional variance explained') + 
      ylab('Percent functional variance explained') +
      facet_wrap(~name, ncol = 5)
  } else if (make_norm & !r2_flag){
    p <- coup_df %>%
      ggplot(aes(x = pc_var, y = value, group = name)) + 
      geom_point() + 
      #geom_path() + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      xlab('Percent compositional variance explained') + 
      ylab('Percent functional variance explained') +
      facet_wrap(~name, ncol = 5)
    
    p2 <- coup_df %>% group_by(pc_var) %>% mutate(mean_value = mean(value)) %>%
      ggplot(aes(x = pc_var, y = mean_value, group = name)) + 
      geom_point() + 
      #geom_path() + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      xlab('Percent compositional variance explained') + 
      ylab('Percent functional variance explained') +
      facet_wrap(~name, ncol = 5)
    
  } else if (r2_flag & !make_norm){
    p <- coup_df %>% 
      ggplot(aes(x = pc_var, y = value*var_res, group = name, color = color)) + 
      geom_point() + 
      #geom_path() + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +  
      scale_color_continuous(type = 'viridis') + 
      xlab('Percent compositional variance explained') + 
      ylab('Percent functional variance explained') +
      facet_wrap(~name, ncol = 5)
  } else { #unnormalized and no color 
    p <- coup_df %>% 
      ggplot(aes(x = pc_var, y = value*var_res, group = name)) + 
      geom_point() + 
      #geom_path() + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab('Percent compositional variance explained') + 
      ylab('Percent functional variance explained') +
      facet_wrap(~name, ncol = 5)
  }
 
  print(p)
  
  return(p)
}

#plot decomp res against thisworks and it works 
thisworks <- coup_full[1,] %*% diag(pc_out$importance[1,])
#sum(thisworks^2) is total variance of a trait 

plot(decomp_res, (thisworks^2)/sum(thisworks^2))
abline(0,1)

#how much variance should pc1 explain
comp_pc1 <- sp_out %*% pc_out$rotation[,1]

summary(lm(Ti_res[,1] ~ comp_pc1))

decomp_res <- c()
for (i in 1:ncol(pc_out$rotation)){
  comp_vec <- sp_out %*% pc_out$rotation[,i]
  decomp_res <- c(decomp_res, summary(lm(Ti_res[,1] ~ comp_vec))$r.squared) 
}

##################
# oof regression #
##################

oof_regression <- function(predictor_res, depend_res, max_npts = 100){
  npts <- length(depend_res)
  
  pred_vals <- c()
  
  if (npts > max_npts){
    npts <- max_npts
  }
  
  for (i in 1:npts){
    y <- as.matrix(depend_res[-i])
    f <- as.formula(y ~ .)
    x <- model.matrix(f, data.frame(predictor_res[-i,]))
    
    topred_f <-  as.formula(as.matrix(1) ~ .)
    topred_x <- model.matrix(topred_f, data.frame(t(as.matrix(predictor_res[i,]))))
    
    cv_fit <- cv.glmnet(x = x, y = y) 
    
    predicted_y <- predict(cv_fit, topred_x, s = "lambda.min")
    pred_vals <- c(pred_vals, predicted_y)
  }
  
  r2 <- get_r2(observed = depend_res[1:npts], predicted = pred_vals) #cor(depend_res, pred_vals)^2
  #plot(depend_res[1:npts], pred_vals)
  #abline(0,1)
  
  return(r2)
}

