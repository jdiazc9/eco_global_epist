rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(cowplot)

# load data sets
files <- list.files('../data_sets', full.names = T)
data <- lapply(files, FUN = function(file) read.csv(file))

# break Kuebbing data into two subsets (invaasive and native)
n_inv <- sapply(2:nrow(data[[4]]),
                FUN = function(i) any(data[[4]][i, 5:8] != 0))

data[[6]] <- data[[4]][c(TRUE, n_inv), c(5:9)]
data[[4]] <- data[[4]][c(TRUE, !n_inv), c(1:4, 9)]

files [6] <- paste(files[4], '_inv', sep = '')

# add our own data (pyoverdine secretion)
data[[7]] <- read.csv('../pyoverdine_data/training_set.csv')
data[[7]]$function. <- rowMeans(data[[7]][, 9:11])
data[[7]] <- data[[7]][, c(1:8, 12)]

test_set <- read.csv('../pyoverdine_data/test_set.csv')
test_set$function. <- rowMeans(test_set[, c(9, 10)])
test_set <- test_set[, c(1:8, 11)]

data[[7]] <- rbind(data[[7]], test_set)

# add genetic data
data[[8]] <- read.csv('../genetic_data_sets/Khan_fitness.csv')

# average function of repeats
data <- lapply(data,
               FUN = function(data_i) aggregate(formula = function. ~., data = data_i, FUN = mean))

# get (empirical) slopes for every species in every data set
ge_data <- lapply(data, FUN = function(x) makeGEdata(matrix2string(x)))
fees <- lapply(ge_data, makeFEEs)
for (i in 1:length(data)) {
  fees[[i]]$species <- rownames(fees[[i]])
  fees[[i]]$dataset <- c(basename(files), 'pyo', 'khan')[i]
}
fees <- do.call(rbind, fees)
fees <- data.frame(dataset = fees$dataset,
                   species = fees$species,
                   slope = fees$b,
                   intercept = fees$a)
fees$dataset[is.na(fees$dataset)] <- 'pyo'

# estimate slopes from pervasive pairwise epistasis
params <- data.frame(dataset = character(0),
                     species_i = character(0),
                     species_j = character(0),
                     deltaF_j_Bi = numeric(0),
                     eps_ij = numeric(0),
                     w_ij = numeric(0),
                     b_ij = numeric(0)) # initialize output data set

for (ds in 1:length(data)) { # loop through data sets
  
  print(paste('DATA SET #', ds, sep = ''))
  
  data_i <- data[[ds]]
  
  sp <- colnames(data_i)[1:(ncol(data_i)-1)] # species namees in dataset
  
  for (i in 1:length(sp)) { # loop through all focal species in dataset
    
    print(paste('   species', i))
    
    sp_i <- sp[i] # species i
    sp_other <- sp[!(sp == sp_i)] # other species
    
    deltaF_j_Bi <- setNames(rep(NA, length(sp_other)), sp_other)
    eps_ij <- setNames(rep(NA, length(sp_other)), sp_other)
    
    for (sp_j in sp_other) { # loop through species j
      
      # deltaF_j_Bi
      B_i <- data_i[data_i[, colnames(data_i) == sp_i] == 0, ] # backgrounds of species i, B(i)
      B_i_j <- B_i[B_i[, colnames(B_i) == sp_j] == 1, ] # backgrounds of i that contain species j
      B_i_noj <- B_i[B_i[, colnames(B_i) == sp_j] == 0, ] # backgrounds of i that do not contain species j
      
      deltaF_j_Bi[sp_j] <- mean(B_i_j$function., na.rm = T) - mean(B_i_noj$function., na.rm = T) # estimated avg. functional effect of species j on the backgrounds of species i
      
      # epsilon_ij
      B_i <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 0 & data_i[colnames(data_i) == sp_j] == 0, ]) # backgrounds of i not containing j
      B_i_j <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 0 & data_i[colnames(data_i) == sp_j] == 1, ]) # backgrounds of i containing j
      B_i_i <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 1 & data_i[colnames(data_i) == sp_j] == 0, ]) # backgrounds of i with no j, with i knocked in
      B_i_ij <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 1 & data_i[colnames(data_i) == sp_j] == 1, ]) # backgrounds of i containing j, with i knocked in
      
      B_i$background <- B_i$community
      B_i_j$background <- sapply(B_i_j$community,
                                 FUN = function(x) {
                                   bg <- strsplit(x, split = ',')[[1]]
                                   bg <- bg[bg != sp_j]
                                   bg <- paste(bg, collapse = ',')
                                   return(bg)
                                 })
      B_i_i$background <- sapply(B_i_i$community,
                                 FUN = function(x) {
                                   bg <- strsplit(x, split = ',')[[1]]
                                   bg <- bg[bg != sp_i]
                                   bg <- paste(bg, collapse = ',')
                                   return(bg)
                                 })
      B_i_ij$background <- sapply(B_i_ij$community,
                                  FUN = function(x) {
                                    bg <- strsplit(x, split = ',')[[1]]
                                    bg <- bg[bg != sp_j & bg != sp_i]
                                    bg <- paste(bg, collapse = ',')
                                    return(bg)
                                  })
      
      allF <- merge(B_i, B_i_i, by = 'background', suffixes = c('_B', '_Bi'))
      allF <- merge(allF, B_i_j, by = 'background')
      colnames(allF)[6:7] <- paste(colnames(allF)[6:7], '_Bj', sep = '')
      allF <- merge(allF, B_i_ij, by = 'background')
      colnames(allF)[8:9] <- paste(colnames(allF)[8:9], '_Bij', sep = '')
      
      allF <- allF[, c('community_B', 'community_Bi', 'community_Bj', 'community_Bij', 'fun_B', 'fun_Bi', 'fun_Bj', 'fun_Bij')]
      
      eps_ij[sp_j] <- mean(allF$fun_Bij - (allF$fun_Bi + allF$fun_Bj - allF$fun_B), na.rm = T) # deviation from additivity
      
    }
    
    #attach to output data set
    params <- rbind(params,
                    data.frame(dataset = c(basename(files), 'pyo', 'khan')[ds],
                               species_i = sp_i,
                               species_j = sp_other,
                               deltaF_j_Bi = as.numeric(deltaF_j_Bi),
                               eps_ij = as.numeric(eps_ij)))
  
  }

}

# P. polymyxa needs a different analysis: backgrounds have to be split between those with/without B. thuringiensis
ds <- 1
data_i <- data[[ds]]
ge_data_i <- makeGEdata(matrix2string(data[[ds]]))
sp <- colnames(data_i)[1:(ncol(data_i)-1)]
sp_i <- 'P'
sp_other <- sp[!(sp %in% c('P', 'T'))]

for (thurin in c(0, 1)) {
  
  # observed slopes & intercepts
  focal_backgrounds <- as.logical(sapply(ge_data_i$background, FUN = function(x) grepl('T', x)))
  if (!thurin) focal_backgrounds <- !focal_backgrounds
  obs_slope <- as.numeric(lm(data = ge_data_i[ge_data_i$knock_in == 'P' & focal_backgrounds, ], formula = d_f ~ background_f)$coefficients[2])
  obs_intercept <- as.numeric(lm(data = ge_data_i[ge_data_i$knock_in == 'P' & focal_backgrounds, ], formula = d_f ~ background_f)$coefficients[1])
  fees <- rbind(fees,
                data.frame(dataset = basename(files)[ds],
                           species = paste('P', thurin, sep = '_'),
                           slope = obs_slope,
                           intercept = obs_intercept))
  
  # parameters to calculate expected slopes
  deltaF_j_Bi <- setNames(rep(NA, length(sp_other)), sp_other)
  eps_ij <- setNames(rep(NA, length(sp_other)), sp_other)
  
  for (sp_j in sp_other) { # loop through species j
    
    # deltaF_j_Bi
    B_i <- data_i[data_i[, colnames(data_i) == sp_i] == 0 & data_i[, colnames(data_i) == 'T'] == thurin, ] # backgrounds of polymyxa that contain or do not contain B. thuringiensis
    B_i_j <- B_i[B_i[, colnames(B_i) == sp_j] == 1, ] # backgrounds of i that contain species j
    B_i_noj <- B_i[B_i[, colnames(B_i) == sp_j] == 0, ] # backgrounds of i that do not contain species j
    
    deltaF_j_Bi[sp_j] <- mean(B_i_j$function., na.rm = T) - mean(B_i_noj$function., na.rm = T) # estimated avg. functional effect of species j on the backgrounds of species i
    
    # epsilon_ij
    B_i <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 0 & data_i[colnames(data_i) == sp_j] == 0 & data_i[colnames(data_i) == 'T'] == thurin, ]) # backgrounds of i not containing j
    B_i_j <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 0 & data_i[colnames(data_i) == sp_j] == 1  & data_i[colnames(data_i) == 'T'] == thurin, ]) # backgrounds of i containing j
    B_i_i <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 1 & data_i[colnames(data_i) == sp_j] == 0  & data_i[colnames(data_i) == 'T'] == thurin, ]) # backgrounds of i with no j, with i knocked in
    B_i_ij <- matrix2string(data_i[data_i[, colnames(data_i) == sp_i] == 1 & data_i[colnames(data_i) == sp_j] == 1  & data_i[colnames(data_i) == 'T'] == thurin, ]) # backgrounds of i containing j, with i knocked in
    
    B_i$background <- B_i$community
    B_i_j$background <- sapply(B_i_j$community,
                               FUN = function(x) {
                                 bg <- strsplit(x, split = ',')[[1]]
                                 bg <- bg[bg != sp_j]
                                 bg <- paste(bg, collapse = ',')
                                 return(bg)
                               })
    B_i_i$background <- sapply(B_i_i$community,
                               FUN = function(x) {
                                 bg <- strsplit(x, split = ',')[[1]]
                                 bg <- bg[bg != sp_i]
                                 bg <- paste(bg, collapse = ',')
                                 return(bg)
                               })
    B_i_ij$background <- sapply(B_i_ij$community,
                                FUN = function(x) {
                                  bg <- strsplit(x, split = ',')[[1]]
                                  bg <- bg[bg != sp_j & bg != sp_i]
                                  bg <- paste(bg, collapse = ',')
                                  return(bg)
                                })
    
    allF <- merge(B_i, B_i_i, by = 'background', suffixes = c('_B', '_Bi'))
    allF <- merge(allF, B_i_j, by = 'background')
    colnames(allF)[6:7] <- paste(colnames(allF)[6:7], '_Bj', sep = '')
    allF <- merge(allF, B_i_ij, by = 'background')
    colnames(allF)[8:9] <- paste(colnames(allF)[8:9], '_Bij', sep = '')
    
    allF <- allF[, c('community_B', 'community_Bi', 'community_Bj', 'community_Bij', 'fun_B', 'fun_Bi', 'fun_Bj', 'fun_Bij')]
    
    eps_ij[sp_j] <- mean(allF$fun_Bij - (allF$fun_Bi + allF$fun_Bj - allF$fun_B), na.rm = T) # deviation from additivity
    
  }
  
  #attach to output data set
  params <- rbind(params,
                  data.frame(dataset = c(basename(files), 'pyo', 'khan')[ds],
                             species_i = paste(sp_i, thurin, sep = '_'),
                             species_j = sp_other,
                             deltaF_j_Bi = as.numeric(deltaF_j_Bi),
                             eps_ij = as.numeric(eps_ij)))
  
}

#params <- params[!(params$dataset == basename(files)[1] & params$species_i == 'P'), ]
#fees <- fees[!(fees$dataset == basename(files)[1] & fees$species == 'P'), ]

# use formula to quantify expected slopes from deltaF and eps_ij measurements
expected_slopes <- data.frame(dataset = character(0),
                              species = character(0),
                              expected_slope = numeric(0))
for (ds in unique(params$dataset)) {
  for (i in unique(params$species_i[params$dataset == ds])) {
    
    params_i <- params[params$dataset == ds & params$species_i == i, ]
    w_ij <- params_i$deltaF_j_Bi^2 / sum(params_i$deltaF_j_Bi^2)
    b_ij <- params_i$eps_ij / params_i$deltaF_j_Bi
    
    expected_slope <- sum(w_ij*b_ij)
    
    if (!is.na(expected_slope)) {
      expected_slopes <- rbind(expected_slopes,
                               data.frame(dataset = ds,
                                          species = i,
                                          expected_slope = expected_slope))
    }
    
  }
}

# use formula to quantify expected intercepts
expected_intercepts <- data.frame(dataset = character(0),
                                  species = character(0),
                                  expected_intercept = numeric(0))
for (ds in 1:length(unique(params$dataset))) {
  
  ge_data_i <- makeGEdata(matrix2string(data[[ds]]))
  
  for (i in unique(params$species_i[params$dataset == unique(params$dataset)[ds]])) {
    
    avg_deltaF_i <- mean(ge_data_i$d_f[ge_data_i$knock_in == i])
    avg_F_Bi <- mean(ge_data_i$background_f[ge_data_i$knock_in == i])
    
    b_i <- expected_slopes$expected_slope[expected_slopes$dataset == unique(params$dataset)[ds] & expected_slopes$species == i]
    
    expected_intercept <- (1 - b_i/2)*avg_deltaF_i - b_i*avg_F_Bi
    
    if (length(expected_intercept)) {
      if (!is.na(expected_intercept)) {
        expected_intercepts <- rbind(expected_intercepts,
                                     data.frame(dataset = unique(params$dataset)[ds],
                                                species = i,
                                                expected_intercept = expected_intercept))
      }
    }
    
  }
}

# again, P. polymyxa requires special treatment (splitting backgrounds by presence/absence of B. thuringiensis)
ds <- 1
data_i <- data[[ds]]
ge_data_i <- makeGEdata(matrix2string(data[[ds]]))

for (thurin in c(0, 1)) {
  
  if (thurin) {
    avg_deltaF_i <- mean(ge_data_i$d_f[ge_data_i$knock_in == 'P' & grepl('T', ge_data_i$background)])
    avg_F_Bi <- mean(ge_data_i$background_f[ge_data_i$knock_in == 'P' & grepl('T', ge_data_i$background)])
  }
  else {
    avg_deltaF_i <- mean(ge_data_i$d_f[ge_data_i$knock_in == 'P' & !grepl('T', ge_data_i$background)])
    avg_F_Bi <- mean(ge_data_i$background_f[ge_data_i$knock_in == 'P' & !grepl('T', ge_data_i$background)])
  }
  
  b_i <- expected_slopes$expected_slope[expected_slopes$dataset == unique(params$dataset)[ds] & expected_slopes$species == paste('P', thurin, sep = '_')]
  
  expected_intercept <- (1 - b_i/2)*avg_deltaF_i - b_i*avg_F_Bi
  
  expected_intercepts <- rbind(expected_intercepts,
                               data.frame(dataset = unique(params$dataset)[ds],
                                          species = paste('P', thurin, sep = '_'),
                                          expected_intercept = expected_intercept))
  
}

# merge data frames and plot
expected_slopes$dataset[expected_slopes$dataset == basename(files)[6]] <- basename(files)[4]
expected_intercepts$dataset[expected_intercepts$dataset == basename(files)[6]] <- basename(files)[4]
fees$dataset[fees$dataset == basename(files)[6]] <- basename(files)[4]


plot_this <- merge(expected_slopes, expected_intercepts, by = c('dataset', 'species'))
plot_this <- merge(fees, plot_this, by = c('dataset', 'species'))

plot_this$dataset <- setNames(c('Bacterial starch hydrolysis',
                                'Bacterial butyrate secretion',
                                'Phytoplankton biomass',
                                'Above-ground plant biomass',
                                'Bacterial xylose oxidation',
                                'Bacterial pyoverdine secretion',
                                'E. coli fitness'),
                              c(basename(files)[1:5], 'pyo', 'khan'))[plot_this$dataset]
plot_this$dataset <- factor(plot_this$dataset, levels = c('E. coli fitness',
                                                          'Above-ground plant biomass',
                                                          'Phytoplankton biomass',
                                                          'Bacterial xylose oxidation',
                                                          'Bacterial starch hydrolysis',
                                                          'Bacterial butyrate secretion',
                                                          'Bacterial pyoverdine secretion'))

plot_this$type <- c('Ecological data set', 'Genetic data set')[1 + plot_this$dataset %in% c('E. coli fitness')]
plot_this$type <- factor(plot_this$type, levels = c('Genetic data set',
                                                    'Ecological data set'))

mycolors <- setNames(c('black',
                       '#d6d62d',
                       '#66b666',
                       '#cb96c3',
                       '#d72027',
                       '#519ed7',
                       'black'),
                     levels(plot_this$dataset))

# standardize intercepts so they can be plotted in a common scale
plot_this$intercept.std <- plot_this$intercept
plot_this$expected_intercept.std <- plot_this$expected_intercept
for (ds in unique(plot_this$dataset)) {
  
  avg <- mean(plot_this$intercept[plot_this$dataset == ds])
  std <- sd(plot_this$intercept[plot_this$dataset == ds])
  
  plot_this$intercept.std[plot_this$dataset == ds] <- (plot_this$intercept[plot_this$dataset == ds] - avg)/std
  plot_this$expected_intercept.std[plot_this$dataset == ds] <- (plot_this$expected_intercept[plot_this$dataset == ds] - avg)/std
  
}

# plot slopes
lims <- c(min(c(plot_this$slope, plot_this$expected_slope)),
          max(c(plot_this$slope, plot_this$expected_slope)))

err <- sd(plot_this$slope - plot_this$expected_slope)

myplot <-
  ggplot(plot_this[!(plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P'), ],
         aes(x = expected_slope, y = slope, color = dataset, shape = type)) +
    geom_abline(slope = 1,
                intercept = 0,
                color = '#d1d3d4') +
    geom_abline(slope = 1,
                intercept = err,
                color = '#d1d3d4',
                linetype = 'dashed') +
    geom_abline(slope = 1,
                intercept = -err,
                color = '#d1d3d4',
                linetype = 'dashed') +
    geom_point(cex = 3) +
    scale_color_manual(values = mycolors,
                       breaks = names(mycolors),
                       name = 'Data set') +
    scale_shape_manual(values = c(16, 3),
                       name = 'Type of data set') +
    scale_x_continuous(name = 'FEE slope\nexpected from pervasive pairwise epistasis',
                       limits = lims) +
    scale_y_continuous(name = 'Observed FEE slope',
                       limits = lims) +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          panel.background = element_blank())

print(myplot)
ggsave(myplot,
       filename = '../plots/slopes_predictions.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# plot intercepts
lims <- c(min(c(plot_this$intercept.std, plot_this$expected_intercept.std)),
          max(c(plot_this$intercept.std, plot_this$expected_intercept.std)))

err <- sd(plot_this$intercept.std - plot_this$expected_intercept.std)

myplot <-
  ggplot(plot_this[!(plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P'), ],
         aes(x = expected_intercept.std, y = intercept.std, color = dataset, shape = type)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_abline(slope = 1,
              intercept = err,
              color = '#d1d3d4',
              linetype = 'dashed') +
  geom_abline(slope = 1,
              intercept = -err,
              color = '#d1d3d4',
              linetype = 'dashed') +
  geom_point(cex = 3) +
  scale_color_manual(values = mycolors,
                     breaks = names(mycolors),
                     name = 'Data set') +
  scale_shape_manual(values = c(16, 3),
                     name = 'Type of data set') +
  scale_x_continuous(name = 'Standardized FEE intercept\nexpected from pervasive pairwise epistasis',
                     limits = lims) +
  scale_y_continuous(name = 'Standardized observed FEE intercept',
                     limits = lims) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank())

print(myplot)
ggsave(myplot,
       filename = '../plots/intercepts_predictions.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# for the xylose dataset, plot FEEs with expected vs. observed fits
ge_data_i <- makeGEdata(matrix2string(data[[5]]))
sp <- unique(ge_data_i$knock_in)

exp_i <- plot_this[plot_this$dataset == 'Bacterial xylose oxidation', c('species', 'expected_slope', 'expected_intercept')]
colnames(exp_i) <- c('knock_in', 'slope', 'intercept')

myplot <-
  ggplot(ge_data_i, aes(x = background_f, y = d_f, group = knock_in)) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_smooth(method = 'lm',
                formula = y~x,
                color = 'black',
                se = T,
                fullrange = T) +
    geom_point(color = 'black',
               cex = 3,
               shape = 16) +
    geom_abline(data = exp_i,
                aes(slope = slope, intercept = intercept),
                linetype = 'dashed',
                color = 'firebrick1',
                size = 1) +
    scale_x_continuous(breaks = pretty_breaks(n = 3),
                       name = 'Function of ecological background [a.u.]') +
    scale_y_continuous(breaks = pretty_breaks(n = 3),
                       name = 'dF [a.u.]') +
    facet_wrap(~knock_in,
               nrow = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

print(myplot)
ggsave(myplot,
       filename = '../plots/xylose_predicted_fees.pdf',
       device = 'pdf',
       dpi = 600,
       width = 250,
       height = 100,
       units = 'mm',
       limitsize = F)

# make supplementary plots for P. polymyxa
ge_data_i <- makeGEdata(matrix2string(data[[1]]))
ge_data_i <- ge_data_i[ge_data_i$knock_in == 'P', ]
ge_data_i$group <- c('A', 'B')[1 + sapply(ge_data_i$background, FUN = function(x) grepl('T', x))]

myplot <-
  ggplot(ge_data_i, aes(x = background_f, y = d_f)) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_smooth(method = 'lm',
                formula = y~x,
                color = 'black',
                se = F,
                fullrange = T) +
    geom_abline(slope = plot_this$expected_slope[plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P'],
                intercept = plot_this$expected_intercept[plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P'],
                linetype = 'dashed',
                color = 'black',
                size = 1) +
    geom_point(color = 'black',
               cex = 3,
               shape = 16) +
    scale_x_continuous(name = 'F (background) [a.u.]',
                       limits = c(-5, 15)) +
    scale_y_continuous(name = 'dF [a.u.]',
                       limits = c(-10, 35)) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

print(myplot)
ggsave(myplot,
       filename = '../plots/polymyxa_predicted_fees_global.pdf',
       device = 'pdf',
       dpi = 600,
       width = 150,
       height = 100,
       units = 'mm',
       limitsize = F)

myplot <-
  ggplot(ge_data_i, aes(x = background_f, y = d_f, group = group, color = group)) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_smooth(method = 'lm',
                formula = y~x,
                se = F,
                fullrange = T) +
    geom_abline(slope = plot_this$expected_slope[plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P_0'],
                intercept = plot_this$expected_intercept[plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P_0'],
                linetype = 'dashed',
                color = '#d72027',
                size = 1) +
    geom_abline(slope = plot_this$expected_slope[plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P_1'],
                intercept = plot_this$expected_intercept[plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P_1'],
                linetype = 'dashed',
                color = '#519ed7',
                size = 1) +
    geom_point(cex = 3,
               shape = 16) +
    scale_x_continuous(name = 'F (background) [a.u.]',
                       limits = c(-5, 15)) +
    scale_y_continuous(name = 'dF [a.u.]',
                       limits = c(-10, 35)) +
    scale_color_manual(values = c('#d72027',
                                  '#519ed7')) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

print(myplot)
ggsave(myplot,
       filename = '../plots/polymyxa_predicted_fees_split.pdf',
       device = 'pdf',
       dpi = 600,
       width = 150,
       height = 100,
       units = 'mm',
       limitsize = F)

  



  