# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(cowplot)
library(ggbreak)

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

for (r in 1:100) { # 1:N randomizations per landscape

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

ggsave(filename = '../plots/benchmark_slope_vs_R2.pdf',
       device = 'pdf',
       dpi = 600,
       width = 330,
       height = 80,
       units = 'mm',
       limitsize = F)

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

ggsave(filename = '../plots/benchmark_slope_vs_intercept.pdf',
       device = 'pdf',
       dpi = 600,
       width = 300,
       height = 80,
       units = 'mm',
       limitsize = F)

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

ggsave(filename = '../plots/benchmark_intercept_vs_R2.pdf',
       device = 'pdf',
       dpi = 600,
       width = 330,
       height = 80,
       units = 'mm',
       limitsize = F)

# statistical tests



