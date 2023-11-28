rm(list = ls())
source('./ecoFunctions.R')
source('./auxFunctions.R')
library(ggh4x)

N <- 6 # number of species
comms <- makeComms(N)

set.seed(1) # for reproducibility

# we model function F(s) of a community s as:
#   F(s) = sum f_i s_i + sum f_ij s_i s_j + alpha EPS(s)
# i.e., additive species contributions, pairwise epistasis, plus a generic "noise"/high-order epistasis term EPS(s) for community s
simulateLandscape <- function(f_i, f_ij, alpha = 0) {
  
  N <- length(f_i)
  comms <- makeComms(N)
  
  f <- sapply(1:nrow(comms),
              FUN = function(comm) {
                
                if(sum(comms[comm, ] == 1) > 1) {
                  
                  pairs <- t(combn(which(comms[comm, ] == 1), 2))
                  return(sum(f_i[comms[comm, ] == 1]) +
                           sum(sapply(1:nrow(pairs),
                                      FUN = function(p) f_ij[pairs[p, 1], pairs[p, 2]])))
                  
                } else {
                  
                  return(sum(f_i[comms[comm, ] == 1]))
                  
                }
                
              })
  
  f <- f + alpha * rnorm(nrow(comms), mean = mean(f), sd = sd(f))
  
  return(cbind(comms, fun = f))
  
}

# make parameters to simulate a collection of landscapes
params_file <- list.files(path = './',
                          pattern = 'simLandscapes_params.txt',
                          full.names = T)

if (length(params_file)) {
  
  params <- read.table(params_file, sep = '\t', header = T)
  
} else {
  
  params <- data.frame(mean_fij = rep(seq(-1, 1, length.out = 7), 7),
                       sd_fij = rep(seq(0.1, 0.9, length.out = 7), each = 7))
  params <- do.call(rbind,
                    lapply(seq(0, 1, length.out = 10),
                           FUN = function(x) cbind(params, alpha = x)))
  params <- do.call(rbind, replicate(5, params, simplify = FALSE))
  
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
              file = './simLandscapes_params.txt',
              quote = F,
              row.names = F,
              col.names = F,
              sep = '\t')
  
  for (p in 1:nrow(params)) {
    
    print(paste('Sim. landscape ', p, ' of ', nrow(params), ' (', round(100*p/nrow(params)), '%)', sep = ''))
    
    # sample coefficients
    f_i <- rnorm(n = N, mean = 0, sd = 1)
    f_ij <- matrix(rnorm(N^2,
                         mean = params$mean_fij[p],
                         sd = params$sd_fij[p]),
                   nrow = N)
    
    # make simulated landscape
    simLandscape <- simulateLandscape(f_i, f_ij, alpha = params$alpha[p])
    
    #plotFitnessGraph(synthLandscape)
    #plotFEEs_clean(synthLandscape)
    
    # get ruggedness
    params$rs[p] <- get_rs(simLandscape)
    
    # evaluate quality of predictions (stitching method) in synthetic landscape
    po <- evaluatePredictions(simLandscape)
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
    po <- get_all_loo_fits(simLandscape)
    
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
                file = './simLandscapes_params.txt',
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t',
                append = T)
    
  }
  
}

### MAKE PLOTS



#heatmaps
mycolors <- c('#99d7dc', '#176766', '#b33a3b')


tst <- aggregate(R2.identity_stitching ~ mean_fij + sd_fij + alpha,
                 data = params,
                 FUN = mean)

ggplot(tst,
       aes(x = mean_fij, y = sd_fij, fill = pmax(R2.identity_stitching, 0))) +
  geom_tile() +
  scale_y_continuous(name = expression(paste('sd ', epsilon[italic(ij)])),
                     expand = c(0, 0),
                     breaks = c(0.1, 0.5, 0.9)) +
  scale_x_continuous(name = expression(paste('mean ', epsilon[italic(ij)])),
                     expand = c(0, 0),
                     breaks = c(-1, 0, 1)) +
  scale_fill_gradient2(name = expression(paste(italic(R)^2, ' predictions vs. observations')),
                       low = '#d32f37',
                       high = '#76d3d6',
                       limits = c(0, 1),
                       midpoint = 0.5,
                       breaks = pretty_breaks(n = 3)) +
  facet_wrap(~ alpha, nrow = 1) +
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
       width = 300,
       height = 125,
       units = 'mm',
       limitsize = F)

# variance fraction explained by HOIs
vH <- do.call(rbind,
              lapply(unique(params$alpha),
                     FUN = function(a) {
                       
                       vH <- sapply(1:100,
                                    FUN = function(i) {
                                      
                                      f_i <- rnorm(n = N, mean = 0, sd = 1)
                                      f_ij <- matrix(rnorm(N^2,
                                                           mean = 0,
                                                           sd = 0.5),
                                                     nrow = N)
                                      
                                      simLandscape <- simulateLandscape(f_i, f_ij, alpha = a)
                                      
                                      return(get_vH(simLandscape))
                                      
                                    })
                       
                       return(data.frame(alpha = a,
                                         vH = vH))
                       
                     }))
vH <- do.call(data.frame,
              aggregate(vH ~ alpha,
                        data = vH,
                        FUN = function(vH) c(mean = mean(vH), sd = sd(vH))))

ggplot(vH, aes(x = alpha, y = 100*vH.mean, ymax = 100*(vH.mean + vH.sd), ymin = 100*(vH.mean - vH.sd))) +
  geom_bar(stat = 'identity',
           fill = 'gray',
           width = 0.05) +
  geom_errorbar(color = 'gray',
                width = 0) +
  scale_y_continuous(name = '% of functional variance\ndue to HOIs',
                     breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(aspect.ratio = 0.1,
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top',
        panel.border = element_blank()) +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/synthLandscapes/vH.pdf',
       device = 'pdf',
       dpi = 600,
       width = 300,
       height = 125,
       units = 'mm',
       limitsize = F)

# same analysis for empirical landscapes

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

emp_vH <- do.call(rbind,
                  lapply(1:length(data),
                         FUN = function(i) {
                           
                           df <- data[[i]]
                           colnames(df)[ncol(df)] <- 'fun'
                           vH <- get_vH(df)
                           
                           return(data.frame(dataset = basename(files)[i],
                                             vH = vH))
                           
                         }))

emp_vH$dataset <- setNames(c('Bacterial starch hydrolysis',
                             'Bacterial butyrate secretion',
                             'Phytoplankton biomass',
                             'Above-ground plant biomass',
                             'Bacterial xylose oxidation',
                             'Bacterial pyoverdine secretion'),
                           basename(files)[1:6])[emp_vH$dataset]
emp_vH$dataset <- factor(emp_vH$dataset, levels = c('Bacterial pyoverdine secretion',
                                                    'Above-ground plant biomass',
                                                    'Phytoplankton biomass',
                                                    'Bacterial xylose oxidation',
                                                    'Bacterial starch hydrolysis',
                                                    'Bacterial butyrate secretion'))

ggplot(emp_vH, aes(x = dataset, y = 100*vH, fill = dataset)) +
  geom_bar(stat = 'identity',
           width = 0.4) +
  scale_y_continuous(name = '% of functional variance\ndue to HOIs',
                     breaks = pretty_breaks(n = 2)) +
  scale_fill_manual(name = '',
                     values = setNames(c('black',
                                         '#d6d62d',
                                         '#66b666',
                                         '#cb96c3',
                                         '#d72027',
                                         '#519ed7'),
                                       levels(emp_vH$dataset))) +
  theme_bw() +
  theme(aspect.ratio = 0.5,
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14,
                                   angle = 30,
                                   hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/synthLandscapes/emp_vH.pdf',
       device = 'pdf',
       dpi = 600,
       width = 150,
       height = 125,
       units = 'mm',
       limitsize = F)

# examine interaction structure (distributions of fitness & epistatic effects) in empirical landscapes

emp_int <- do.call(rbind,
                  lapply(1:length(data),
                         FUN = function(i) {
                           
                           df <- data[[i]]
                           colnames(df)[ncol(df)] <- 'fun'
                           dee <- getDEE(df)
                           dee <- dee[!is.na(dee$eps_ij), ]
                           
                           return(cbind(dataset = basename(files)[i],
                                        dee))
                           
                         }))

emp_int$dataset <- setNames(c('Bacterial starch hydrolysis',
                             'Bacterial butyrate secretion',
                             'Phytoplankton biomass',
                             'Above-ground plant biomass',
                             'Bacterial xylose oxidation',
                             'Bacterial pyoverdine secretion'),
                           basename(files)[1:6])[emp_int$dataset]
emp_int$dataset <- factor(emp_int$dataset, levels = c('Bacterial pyoverdine secretion',
                                                    'Above-ground plant biomass',
                                                    'Phytoplankton biomass',
                                                    'Bacterial xylose oxidation',
                                                    'Bacterial starch hydrolysis',
                                                    'Bacterial butyrate secretion'))

# full species names
sp_names <- vector(mode = 'list', length = 6)
sp_names[[1]] <- setNames(c('B. cereus', 'B. megaterium', 'B. mojavensis', 'P. polymyxa', 'B. subtilis', 'B. thuringiensis'),
                          c('C', 'E', 'M', 'P', 'S', 'T'))
sp_names[[2]] <- setNames(c('P. copri','P. johnsonii','B. vulgatus','B. fragilis','B. ovatus','B. thetaiotaomicron','B. caccae','B. cellulosilyticus','B. uniformis','D. piger','B. longum','B. adolescentis','B. pseudocatenulatum','C. aerofaciens','E. lenta','F. prausnitzii','C. hiranonis','A. caccae','B. hydrogenotrophica','C. asparagiforme','E. rectale','R. intestinalis','C. comes','D. longicatena','D. formicigenerans'),
                          c('PC','PJ','BV','BF','BO','BT','BC','BY','BU','DP','BL','BA','BP','CA','EL','FP','CH','AC','BH','CG','ER','RI','CC','DL','DF'))
sp_names[[3]] <- setNames(c('A. carterae','Tetraselmis sp.','D. tertiolecta','Synechococcus sp.','T. lutea'),
                          c('A','T','D','S','Ti'))
sp_names[[4]] <- setNames(c('A. millefolium','L. capitata','P. virginianum','S. nutans','L. vulgare','L. cuneata','P. vulgaris','P. pratense'),
                          c('as.nat','fa.nat','la.nat','po.nat','as.inv','fa.inv','la.inv','po.inv'))
sp_names[[5]] <- setNames(c('Rhodoferax sp.','Flavobacterium sp.','Sphingoterrabacterium sp.','Burkholderia sp.','S. yanoikuyae','Bacteroidetes sp.'),
                          c('SL68','SL104','SL106','SL187','SL197','SLWC2'))
sp_names[[6]] <- setNames(c('Enterobacter sp.', 'Pseudomonas sp. 02', 'Klebsiella sp.', 'Pseudomonas sp. 03', 'Pseudomonas sp. 04', 'Raoultella sp.', 'Pseudomonas sp. 01', 'Pseudomonas sp. 05'),
                          paste('sp', 1:8, sep = '_'))

for (ds in 1:6) {
  emp_int$species_i[emp_int$dataset == unique(emp_int$dataset)[ds]] <- sp_names[[ds]][emp_int$species_i[emp_int$dataset == unique(emp_int$dataset)[ds]]]
  emp_int$species_j[emp_int$dataset == unique(emp_int$dataset)[ds]] <- sp_names[[ds]][emp_int$species_j[emp_int$dataset == unique(emp_int$dataset)[ds]]]
}




for (ds in 1:length(unique(emp_int$dataset))) {
  
  # distribution of epistatic effects
  tst <- emp_int[emp_int$dataset == unique(emp_int$dataset)[ds], ]
  indx_i <- sapply(tst$species_i,
                   FUN = function(sp) which(unique(tst$species_i) == sp))
  indx_j <- sapply(tst$species_j,
                   FUN = function(sp) which(unique(tst$species_i) == sp))
  tst <- tst[indx_j > indx_i, ]
  tst$species_i <- factor(tst$species_i,
                          levels = unique(c(tst$species_i, tst$species_j)))
  tst$species_j <- factor(tst$species_j,
                          levels = unique(c(tst$species_i, tst$species_j)))
  
  ggplot(tst, aes(x = 0, y = eps_ij)) +
    geom_hline(yintercept = 0,
               color = 'gray') +
    geom_violin(fill = 'black',
                color = NA,
                alpha = 0.25,
                width = 0.6,
                scale = 'width') +
    geom_point() +
    facet_grid2(species_i ~ species_j,
                render_empty = F) +
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    scale_y_continuous(name = expression(paste('Interaction between species ', italic(i), ' and ', italic(j), ', ', italic(epsilon)[italic(ij)])),
                       breaks = pretty_breaks(n = 3)) +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = 'none',
          strip.background = element_blank(),
          strip.text = element_text(size = 12,
                                    face = 'italic'),
          strip.text.x = element_text(angle = 90,
                                      hjust = 0),
          strip.text.y = element_text(angle = 0,
                                      hjust = 0),
          strip.clip = 'off')
  
  b <- (400 - 125) / (25 - 6)
  a <- 125 - 6*b + 50
    
  ggsave(filename = paste('../plots/synthLandscapes/distEpistEffects_', basename(files)[ds], '.pdf', sep = ''),
         device = 'pdf',
         dpi = 600,
         width = a + b * length(unique(c(tst$species_i, tst$species_j))),
         height = a + b * length(unique(c(tst$species_i, tst$species_j))),
         units = 'mm',
         limitsize = F)
  
  
  # dsitribution of functional effects
  dfe <- makeGEdata(matrix2string(data[[ds]]))
  dfe$knock_in <- sp_names[[ds]][dfe$knock_in]
  
  ggplot(dfe, aes(x = knock_in, y = d_f)) +
    geom_hline(yintercept = 0,
               color = 'gray') +
    geom_violin(fill = 'black',
                color = NA,
                alpha = 0.25,
                width = 0.6,
                scale = 'width') +
    geom_point() +
    scale_x_discrete() +
    scale_y_continuous(name = expression(paste('Functional effect, ', Delta*italic(F))),
                       breaks = pretty_breaks(n = 3)) +
    theme_bw() +
    theme(aspect.ratio = 0.2,
          panel.grid = element_blank(),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 30,
                                     size = 12,
                                     hjust = 1),
          axis.title.x = element_blank(),
          legend.position = 'none')
  
  ggsave(filename = paste('../plots/synthLandscapes/distFunEffects_', basename(files)[ds], '.pdf', sep = ''),
         device = 'pdf',
         dpi = 600,
         width = 280,
         height = 300,
         units = 'mm',
         limitsize = F)
  
}





