# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(cowplot)
library(ggbreak)
library(tidyr)
library(plotly)

saveplots <- F # set to TRUE to save generated plots under ./plots directory

# load data sets
files <- list.files('../data_sets', full.names = T)
files <- c(files, '../pyoverdine_data/training_set.csv')
#files <- c(files, '../genetic_data_sets/Khan_fitness.csv')
data <- lapply(files, FUN = function(file) read.csv(file))
data[[3]]$function. <- data[[3]]$function. / 1e4
data[[6]]$function. <- rowMeans(data[[6]][, 9:11])
data[[6]] <- data[[6]][, c(1:8, 12)]
data <- lapply(data, FUN = function(df) aggregate(fun ~ ., data = matrix2string(df), FUN = mean))
# data <- do.call(rbind,
#                 lapply(1:length(data), FUN = function(i) cbind(dataset = basename(files)[i], data[[i]])))

# full species names
sp_names <- vector(mode = 'list', length = length(data))

sp_names[[1]] <- setNames(c('B. cereus', 'B. megaterium', 'B. mojavensis', 'P. polymyxa', 'B. subtilis', 'B. thuringiensis'),
                          c('C', 'E', 'M', 'P', 'S', 'T'))
sp_names[[2]] <- setNames(c('P. copri','P. johnsonii','B. vulgatus','B. fragilis','B. ovatus','B. thetaiotaomicron','B. caccae','B. cellulosilyticus','B. uniformis','D. piger','B. longum','B. adolescentis','B. pseudocatenulatum','C. aerofaciens','E. lenta','F. prausnitzii','C. hiranonis','A. caccae','B. hydrogenotrophica','C. asparagiforme','E. rectale','R. intestinalis','C. comes','D. longicatena','D. formicigenerans'),
                          c('PC','PJ','BV','BF','BO','BT','BC','BY','BU','DP','BL','BA','BP','CA','EL','FP','CH','AC','BH','CG','ER','RI','CC','DL','DF'))
sp_names[[3]] <- setNames(c('A. carterae','Tetraselmis sp.','D. tertiolecta','Synechococcus sp.','T. lutea'),
                          c('A','T','D','S','Ti'))
sp_names[[4]] <- setNames(c('A. millefolium','L. capitata','P. virginianum','S. nutans','L. vulgare','L. cuneata','P. vulgaris','P.pratense'),
                          c('as.nat','fa.nat','la.nat','po.nat','as.inv','fa.inv','la.inv','po.inv'))
sp_names[[5]] <- setNames(c('Rhodoferax sp.','Flavobacterium sp.','Sphingoterrabacterium sp.','Burkholderia sp.','S. yanoikuyae','Bacteroidetes sp.'),
                          c('SL68','SL104','SL106','SL187','SL197','SLWC2'))
sp_names[[6]] <- setNames(c('Enterobacter sp.', 'Pseudomonas sp. 02', 'Klebsiella sp.', 'Pseudomonas sp. 03', 'Pseudomonas sp. 04', 'Raoultella sp.', 'Pseudomonas sp. 01', 'Pseudomonas sp. 05'),
                          paste('sp', 1:8, sep = '_'))
sp_names[[7]] <- setNames(c('glmUS', 'pykF', 'rbs', 'spoT', 'topA'),
                          c('g', 'p', 'r', 's', 't'))

randomizeLandscape <- function(df) {
  df$fun <- sample(df$fun)
  return(df)
}

additiveLandscape <- function(df) {
  
  gedf <- makeGEdata(df)
  f_add <- aggregate(d_f ~ knock_in,
                     data = gedf,
                     FUN = mean)
  f_add <- setNames(f_add$d_f, f_add$knock_in)
  additive_model <- sapply(df$community,
                           FUN = function(comm) {
                             comm <- strsplit(comm, split = ',')[[1]]
                             if (length(comm)) {
                               fun <- sum(f_add[comm])
                             } else {
                               fun <- 0
                             }
                             return(fun)
                           })
  return(data.frame(community = df$community,
                    fun = additive_model))
  
}

# empirical FEEs
empirical_fees <- do.call(rbind,
                          lapply(1:length(data),
                                 FUN = function(i) {
                                   gedf <- makeGEdata(data[[i]])
                                   fees <- makeFEEs(gedf)
                                   R2 <- sapply(rownames(fees),
                                                FUN = function(sp) {
                                                  lmod <- lm(d_f ~ background_f,
                                                             data = gedf[gedf$knock_in == sp, ])
                                                  return(summary(lmod)$r.squared)
                                                })
                                   return(data.frame(alpha = NA,
                                                     dataset = basename(files)[i],
                                                     species = rownames(fees),
                                                     slope = fees$b,
                                                     intercept = fees$a,
                                                     R2 = as.numeric(R2)))
                                 }))
     

# go through datasets and randomize
randomizations <- data.frame(rnd_id = numeric(0),
                             alpha = numeric(0),
                             dataset = character(0),
                             species = character(0),
                             slope = numeric(0),
                             intercept = numeric(0),
                             R2 = numeric(0))

for (r in 1:1000) { # 1:N randomizations per landscape

  for (i in 1:length(data)) {
    
    print(paste('Randomization #', r, ', dataset #', i, sep = ''))
  
    df <- data[[i]]
    df_rnd <- randomizeLandscape(df)
    df_add <- additiveLandscape(df)
    
    alpha <- runif(1)# alpha <- sample(c(0, 1), size = 1) # alpha <- 0 #  # this parameter controls whether we benchmark against the additive model (alpha = 1) or the randomized model (alpha = 0)
    df_benchmark <- data.frame(community = df$community,
                               fun = alpha*df_add$fun + (1-alpha)*df_rnd$fun)
    
    gedf <- makeGEdata(df_benchmark)
    fees_bench <- makeFEEs(gedf)
    R2 <- sapply(rownames(fees_bench),
                 FUN = function(sp) {
                   lmod <- lm(d_f ~ background_f,
                              data = gedf[gedf$knock_in == sp, ])
                   return(summary(lmod)$r.squared)
                 })
    
    randomizations <- rbind(randomizations,
                            data.frame(rnd_id = r,
                                       alpha = alpha,
                                       dataset = basename(files)[i],
                                       species = rownames(fees_bench),
                                       slope = fees_bench$b,
                                       intercept = fees_bench$a,
                                       R2 = as.numeric(R2)))
  
  }
}

# plot
randomizations$dataset <- factor(randomizations$dataset,
                                 levels = c("training_set.csv",
                                            "plant-biommass_Kuebbing2016_all.csv",
                                            "phytoplankton-biomass_Ghedini2022.csv",
                                            "xylose_Langenheder2010.csv",
                                            "amyl_Sanchez-Gorostiaga2019.csv",
                                            "butyrate_Clark2021.csv"))
empirical_fees$dataset <- factor(empirical_fees$dataset,
                                 levels = levels(randomizations$dataset))

ggplot(randomizations, aes(x = slope, y = R2, color = alpha)) +
  geom_point(shape = 19,
             alpha = 0.5,
             stroke = 0,
             cex = 1.5) +
  geom_point(data = empirical_fees,
             color = 'black',
             shape = 15) +
  facet_wrap(~ dataset,
             nrow = 1,
             scales = 'free') +
  theme_bw() +
  scale_x_continuous(name = 'FEE slope',
                     breaks = pretty_breaks(n = 2),
                     limits = c(-2.5, 1)) +
  scale_y_continuous(name = expression(paste(italic(R)^2, ' of FEE', sep = '')),
                     breaks = pretty_breaks(n = 2),
                     limits = c(0, 1)) +
  scale_color_gradient(low = '#d32f37',
                       high = '#76d3d6') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.background = element_blank())

if(saveplots) {
  ggsave(filename = '../plots/benchmark_slope_vs_R2.pdf',
       device = 'pdf',
       dpi = 600,
       width = 330,
       height = 80,
       units = 'mm',
       limitsize = F)
}


ggplot(randomizations, aes(x = slope, y = intercept, color = alpha)) +
  geom_point(shape = 19,
             alpha = 0.5,
             stroke = 0,
             cex = 1.5) +
  geom_point(data = empirical_fees,
             color = 'black',
             shape = 15) +
  facet_wrap(~ dataset,
             nrow = 1,
             scales = 'free') +
  theme_bw() +
  scale_x_continuous(name = 'FEE slope',
                     breaks = pretty_breaks(n = 2),
                     limits = c(-2.5, 1)) +
  scale_y_continuous(name = 'FEE intercept [a.u.]',
                     breaks = pretty_breaks(n = 2)) +
  scale_color_gradient(low = '#d32f37',
                       high = '#76d3d6') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.background = element_blank())

if(saveplots) {
  ggsave(filename = '../plots/benchmark_slope_vs_intercept.pdf',
         device = 'pdf',
         dpi = 600,
         width = 300,
         height = 80,
         units = 'mm',
         limitsize = F)
}

ggplot(randomizations, aes(x = intercept, y = R2, color = alpha)) +
  geom_point(shape = 19,
             alpha = 0.5,
             stroke = 0,
             cex = 1.5) +
  geom_point(data = empirical_fees,
             color = 'black',
             shape = 15) +
  facet_wrap(~ dataset,
             nrow = 1,
             scales = 'free') +
  theme_bw() +
  scale_x_continuous(name = 'FEE intercept [a.u.]',
                     breaks = pretty_breaks(n = 2)) +
  scale_y_continuous(name = expression(paste(italic(R)^2, ' of FEE', sep = '')),
                     breaks = pretty_breaks(n = 2),
                     limits = c(0, 1)) +
  scale_color_gradient(low = '#d32f37',
                       high = '#76d3d6') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.background = element_blank())

if(saveplots) {
  ggsave(filename = '../plots/benchmark_intercept_vs_R2.pdf',
       device = 'pdf',
       dpi = 600,
       width = 330,
       height = 80,
       units = 'mm',
       limitsize = F)
}




### PROJECTION PLOTS

# wrapper function to make plot from a data frame
makeProjPlot <- function(df, showaxes = F) {
  
  colnames(df) <- c('x', 'y', 'z', 'color')
  
  mins <- apply(df[, 1:3], FUN = min, MARGIN = 2)
  maxs <- apply(df[, 1:3], FUN = max, MARGIN = 2)
  
  expands <- apply(df[, 1:3],
                   FUN = function(x) diff(range(x))*0.1,
                   MARGIN = 2)
  mins <- mins - expands
  maxs <- maxs + expands
  
  pplot <- plot_ly(x = ~x, y = ~y, z = ~z,
                   marker = list(size = 6,
                                 line = list(width = 0),
                                 colorscale = list(c(0, 1), c('#d32f37', '#76d3d6')))) %>%
    layout(scene = list(camera = list(projection = list(type = 'orthographic')),
                        aspectmode = 'manual',
                        aspectratio = list(x = 1, y = 1, z = 1),
                        xaxis = list(range = c(mins['x'], maxs['x']),
                                     title = c('', 'x')[1 + showaxes],
                                     zeroline = F,
                                     showticklabels = showaxes),
                        yaxis = list(range = c(mins['y'], maxs['y']),
                                     title = c('', 'y')[1 + showaxes],
                                     zeroline = F,
                                     showticklabels = showaxes),
                        zaxis = list(range = c(mins['z'], maxs['z']),
                                     title = c('', 'z')[1 + showaxes],
                                     zeroline = F,
                                     showticklabels = showaxes)),
           showlegend = F)
  
  for (i in 1:3) {
    
    dfi <- df
    dfi[, i] <- mins[i]
    
    pplot <- pplot %>% add_markers(data = dfi[1:(nrow(dfi) - 1), ],
                                   marker = list(color = ~color),
                                   opacity = 0.66)
    pplot <- pplot %>% hide_colorbar()
    pplot <- pplot %>% add_markers(data = dfi[nrow(dfi), ],
                                   opacity = 1,
                                   marker = list(color = 'black',
                                                 size = 10,
                                                 symbol = 'square'))
    pplot <- pplot %>% hide_colorbar()
    
  }
  
  manual_axes <- data.frame(x = c(mins['x'], maxs['x'], rep(mins['x'], 4)),
                            y = c(rep(mins['y'], 3), maxs['y'], rep(mins['y'], 2)),
                            z = c(rep(mins['z'], 5), maxs['z']))
  
  pplot <- pplot %>% add_paths(data = manual_axes, x = ~x, y = ~y, z = ~z,
                               marker = list(size = 5,
                                             color = 'black',
                                             opacity = 0.01),
                               line = list(width = 4,
                                           color = 'black'))
  
  return(pplot)
  
}

for (i in 1:length(data)) {
  
  print(paste('DATASET #', i, sep = ''))
  dfi <- rbind(randomizations[randomizations$dataset == basename(files)[i], c(2:ncol(randomizations))],
               empirical_fees[empirical_fees$dataset == basename(files)[i], ])
  
  for (sp in unique(dfi$species)) {
    
    print(paste('   species #', which(sp == unique(dfi$species)), sep = ''))
    
    dfisp <- dfi[dfi$species == sp, c('slope', 'intercept', 'R2', 'alpha')]
    
    pplot <- makeProjPlot(dfisp)
    if(saveplots) {
      orca(pplot,
           file = paste('../plots/projections/',
                        gsub('\\.csv', '', basename(files)[i]),
                        '/',
                        sp_names[[i]][sp],
                        '.pdf', sep = ''),
           scale = 2)
    }
    
  }
  
}





### STATISTICAL TESTS
# mytests <- do.call(rbind,
#                    lapply(1:length(files),
#                           FUN = function(i) {
#                             
#                             out <- do.call(rbind,
#                                            lapply(unique(randomizations$species[randomizations$dataset == basename(files)[i]]),
#                                                   FUN = function(sp) {
#                                                     
#                                                     df <- rbind(randomizations[randomizations$dataset == basename(files)[i] & randomizations$species == sp, c('slope', 'intercept', 'R2')],
#                                                                 empirical_fees[empirical_fees$dataset == basename(files)[i] & empirical_fees$species == sp, c('slope', 'intercept', 'R2')])
#                                                     
#                                                     mytest1 <- ks.test(df$slope[1:(nrow(df)-1)],
#                                                                        df$slope[nrow(df)])
#                                                     mytest2 <- ks.test(df$intercept[1:(nrow(df)-1)],
#                                                                        df$intercept[nrow(df)])
#                                                     mytest3 <- ks.test(df$R2[1:(nrow(df)-1)],
#                                                                        df$R2[nrow(df)])
#                                                     
#                                                     return(data.frame(dataset = basename(files)[i],
#                                                                       species = sp,
#                                                                       pval_slope = mytest1$p.value,
#                                                                       pval_intercept = mytest2$p.value,
#                                                                       pval_R2 = mytest3$p.value))
#                                                     
#                                                   }))
#                             
#                             return(out)
#                             
#                           }))
# 
# pt <- 0.01
# sum(mytests$pval_slope < pt | mytests$pval_intercept < pt | mytests$pval_R2 < pt) / nrow(mytests)

mytests <- do.call(rbind,
                   lapply(1:length(files),
                          FUN = function(i) {
                            
                            out <- do.call(rbind,
                                           lapply(unique(randomizations$species[randomizations$dataset == basename(files)[i]]),
                                                  FUN = function(sp) {
                                                    
                                                    df <- rbind(randomizations[randomizations$dataset == basename(files)[i] & randomizations$species == sp, c('slope', 'intercept', 'R2')],
                                                                empirical_fees[empirical_fees$dataset == basename(files)[i] & empirical_fees$species == sp, c('slope', 'intercept', 'R2')])
                                                    df <- as.data.frame(scale(df))
                                                    
                                                    d <- 0.1 # define interval for choosing null models to compare to
                                                    
                                                    # slope test
                                                    dfi <- df[df$intercept > df$intercept[nrow(df)] - d & df$intercept < df$intercept[nrow(df)] + d &
                                                                df$R2 > df$R2[nrow(df)] - d & df$R2 < df$R2[nrow(df)] + d, ]
                                                    if(nrow(dfi) == 1) {
                                                      pval_slope <- 1/nrow(df) # if no null models yield comparable values, the probability of the empirical FEE resulting from a null model is lower than 1/N (we set it to 1/N)
                                                    } else {
                                                      mytest <- ks.test(dfi$slope[1:(nrow(dfi) - 1)], y = dfi$slope[nrow(dfi)])
                                                      pval_slope <- mytest$p.value
                                                    }
                                                    
                                                    # intercept test
                                                    dfi <- df[df$slope > df$slope[nrow(df)] - d & df$slope < df$slope[nrow(df)] + d &
                                                                df$R2 > df$R2[nrow(df)] - d & df$R2 < df$R2[nrow(df)] + d, ]
                                                    if(nrow(dfi) == 1) {
                                                      pval_intercept <- 1/nrow(df) # same as before
                                                    } else {
                                                      mytest <- ks.test(dfi$intercept[1:(nrow(dfi) - 1)], y = dfi$intercept[nrow(dfi)])
                                                      pval_intercept <- mytest$p.value
                                                    }
                                                    
                                                    # R2 test
                                                    dfi <- df[df$slope > df$slope[nrow(df)] - d & df$slope < df$slope[nrow(df)] + d &
                                                                df$intercept > df$intercept[nrow(df)] - d & df$intercept < df$intercept[nrow(df)] + d, ]
                                                    if(nrow(dfi) == 1) {
                                                      pval_R2 <- 1/nrow(df) # same as before
                                                    } else {
                                                      mytest <- ks.test(dfi$R2[1:(nrow(dfi) - 1)], y = dfi$R2[nrow(dfi)])
                                                      pval_R2 <- mytest$p.value
                                                    }
                                                    
                                                    return(data.frame(dataset = basename(files)[i],
                                                                      species = sp,
                                                                      pval_slope = pval_slope,
                                                                      pval_intercept = pval_intercept,
                                                                      pval_R2 = pval_R2))
                                                    
                                                  }))
                            
                            return(out)
                            
                          }))

pt <- 0.01
sum(mytests$pval_slope < pt | mytests$pval_intercept < pt | mytests$pval_R2 < pt) / nrow(mytests)








