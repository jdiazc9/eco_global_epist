# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(cowplot)
library(tidyr)
library(CEGO)

# load data sets
files <- list.files('../data_sets', full.names = T)
data <- lapply(files, FUN = function(file) read.csv(file))

# break Kuebbing data into two subsets (invasive and native)
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

names(data) <- c(basename(files), 'pyo', 'khan')

# average function of repeats
data <- lapply(data,
               FUN = function(data_i) aggregate(function. ~., data = data_i, FUN = mean))

# scale values for phytoplankton data
data[[3]]$function. <- data[[3]]$function./1e4

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
fees$dataset <- gsub('_inv', '', fees$dataset)

######################################################
### estimate slopes from pervasive pairwise epistasis
######################################################

params <- data.frame(dataset = character(0),
                     background = character(0),
                     species_i = character(0),
                     species_j = character(0),
                     deltaF_j = numeric(0),
                     eps_ij = numeric(0),
                     w_ij = numeric(0),
                     b_ij = numeric(0)) # initialize output data set

for (ds in 1:length(data)) { # loop through data sets
  
  print(paste('DATA SET #', ds, sep = ''))
  
  data_i <- data[[ds]]
  
  sp <- colnames(data_i)[1:(ncol(data_i)-1)] # species namees in dataset
  
  F0 <- data_i$function.[rowSums(data_i[, 1:(ncol(data_i) - 1)]) == 0] # function of 'empty' community
  
  for (i in 1:length(sp)) { # loop through all focal species in dataset
    
    print(paste('   species', i))
    
    sp_i <- sp[i] # species i
    sp_other <- sp[!(sp == sp_i)] # other species
    
    for (sp_j in sp_other) { # loop through species j
      
      # deltaF_j
      deltaF_j <- data_i$function.[rowSums(data_i[, 1:(ncol(data_i) - 1)]) == 1 & data_i[, colnames(data_i) == sp_j] == 1] - F0 # functional effect of species j on 'empty' background
      
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
      allF <- allF[allF$community_B == '', ]
      
      eps_ij <- allF$fun_Bij - (allF$fun_Bi + allF$fun_Bj - allF$fun_B) # deviation from additivity
      if (!length(eps_ij)) eps_ij <- NA
      if (!length(deltaF_j)) deltaF_j <- NA
      
      #attach to output data set
      params <- rbind(params,
                      data.frame(dataset = c(basename(files), 'pyo', 'khan')[ds],
                                 background = '',
                                 species_i = sp_i,
                                 species_j = sp_j,
                                 deltaF_j = as.numeric(deltaF_j),
                                 eps_ij = as.numeric(eps_ij)))
      
    }
  
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
  F0 <- data_i$function.[rowSums(data_i[, 1:(ncol(data_i) - 1)]) == thurin & data_i[, colnames(data_i) == 'T'] == thurin] # functional effect of species j on 'empty' background
  
  for (sp_j in sp_other) { # loop through species j
    
    # deltaF_j
    deltaF_j <- data_i$function.[rowSums(data_i[, 1:(ncol(data_i) - 1)]) == (1 + thurin) & data_i[, colnames(data_i) == sp_j] == 1 & data_i[, colnames(data_i) == 'T'] == thurin] - F0 # functional effect of species j on 'empty' background
    
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
    allF <- allF[allF$community_B == c('', 'T')[1 + thurin], ]
    
    eps_ij <- allF$fun_Bij - (allF$fun_Bi + allF$fun_Bj - allF$fun_B) # deviation from additivity
    
    if (!length(eps_ij)) eps_ij <- NA
    if (!length(deltaF_j)) deltaF_j <- NA
    
    #attach to output data set
    params <- rbind(params,
                    data.frame(dataset = c(basename(files), 'pyo', 'khan')[ds],
                               background = c('', 'T')[1 + thurin],
                               species_i = paste(sp_i, thurin, sep = '_'),
                               species_j = sp_j,
                               deltaF_j = as.numeric(deltaF_j),
                               eps_ij = as.numeric(eps_ij)))
    
  }
  
}

# plot effect sizes
plot_this <- params[params$dataset != 'khan', ]
plot_this$dataset <- gsub('_inv', '', plot_this$dataset)

plot_this$dataset <- setNames(c('Bacterial starch hydrolysis',
                                'Bacterial butyrate secretion',
                                'Phytoplankton biomass',
                                'Above-ground plant biomass',
                                'Bacterial xylose oxidation',
                                'Bacterial pyoverdine secretion'),
                              c(basename(files)[1:5], 'pyo'))[plot_this$dataset]
plot_this$dataset <- factor(plot_this$dataset, levels = c('Above-ground plant biomass',
                                                          'Phytoplankton biomass',
                                                          'Bacterial xylose oxidation',
                                                          'Bacterial starch hydrolysis',
                                                          'Bacterial butyrate secretion',
                                                          'Bacterial pyoverdine secretion'))

sp_names <- vector(mode = 'list', length = 6)
sp_names[[1]] <- setNames(c('B. cereus', 'B. megaterium', 'B. mojavensis', 'P. polymyxa', 'P. polymyxa\n(B. thuringiensis\nnot in background)', 'P. polymyxa\n(B. thuringiensis\nin background)', 'B. subtilis', 'B. thuringiensis'),
                          c('C', 'E', 'M', 'P', 'P_0', 'P_1', 'S', 'T'))
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
sp_names <- unlist(sp_names)

sp_names_abbr <- gsub('\\.[[:space:]]', '', sp_names)
sp_names_abbr <- substr(sp_names_abbr, 1, 2)
sp_names_abbr[2] <- 'Bme'
sp_names_abbr[3] <- 'Bmo'
sp_names_abbr[15] <- 'Bca'
sp_names_abbr[16] <- 'Bce'
sp_names_abbr[22] <- 'Cae'
sp_names_abbr[28] <- 'Cas'
sp_names_abbr[26] <- 'Acac'
sp_names_abbr[34] <- 'Acar'
sp_names_abbr[40] <- 'Lca'
sp_names_abbr[44] <- 'Lcu'
sp_names_abbr[41] <- 'Pvi'
sp_names_abbr[45] <- 'Pvu'
sp_names_abbr[c(54, 56, 57, 59, 60)] <- paste('Ps', c('02', '03', '04', '01', '05'), sep = '')


plot_this$species_i <- sp_names[plot_this$species_i]
plot_this$species_j <- sp_names_abbr[plot_this$species_j]

plot_this$product <- abs(plot_this$deltaF_j * plot_this$eps_ij)
plot_this$product_sign <- as.character(sign(plot_this$deltaF_j * plot_this$eps_ij))

plot_this <- plot_this[!is.na(plot_this$product), ]

# plot starch dataset (for main fig)
ggplot(plot_this[plot_this$dataset == 'Bacterial starch hydrolysis', ],
       aes(x = deltaF_j, y = eps_ij, size = product, color = product_sign)) +
  geom_vline(xintercept = 0,
             linetype = 'dashed') +
  geom_hline(yintercept = 0,
             linetype = 'dashed') +
  # geom_point(alpha = 0.5,
  #            shape = 16) +
  geom_text(aes(label = species_j),
            fontface = 'bold.italic') +
  facet_wrap(~ species_i,
             ncol = 6) +
  scale_x_continuous(name = expression(paste(Delta, italic(F)[italic(j)]^0, sep = '')),
                     breaks = pretty_breaks(n = 2),
                     expand = c(0.3, 0.3)) +
  scale_y_continuous(name = expression(paste(italic(epsilon)[italic(ij)]^0, sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = c(0.3, 0.3)) +
  scale_size(range = c(3, 8)) +
  scale_color_manual(values = setNames(c('firebrick1', 'black', 'deepskyblue'),
                                       c('-1', '0', '1'))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 12,
                                  vjust = 0,
                                  angle = 0),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = 'none')

ggsave(filename = '../plots/interaction_map_starch_noavg.pdf',
       device = 'pdf',
       dpi = 600,
       width = 350,
       height = 800,
       units = 'mm',
       limitsize = F)

# plot all datasets (for supplementary)
datasets <- unique(plot_this$dataset)

for (ds in datasets) {
  
  ggplot(plot_this[plot_this$dataset == ds, ],
         aes(x = deltaF_j, y = eps_ij, size = product, color = product_sign)) +
    geom_vline(xintercept = 0,
               linetype = 'dashed') +
    geom_hline(yintercept = 0,
               linetype = 'dashed') +
    # geom_point(alpha = 0.5,
    #            shape = 16) +
    geom_text(aes(label = species_j),
              fontface = 'bold.italic') +
    facet_wrap(~ species_i,
               ncol = 6) +
    scale_x_continuous(name = expression(paste(Delta, italic(F)[italic(j)]^0, sep = '')),
                       breaks = pretty_breaks(n = 2),
                       expand = c(0.3, 0.3)) +
    scale_y_continuous(name = expression(paste(italic(epsilon)[italic(ij)]^0, sep = '')),
                       breaks = pretty_breaks(n = 3),
                       expand = c(0.3, 0.3)) +
    scale_size(range = c(3, 8)) +
    scale_color_manual(values = setNames(c('firebrick1', 'black', 'deepskyblue'),
                                         c('-1', '0', '1'))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 12,
                                    vjust = 0,
                                    angle = 0),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          legend.position = 'none')
  
  ggsave(filename = paste('../plots/interaction_map_', gsub(' ', '', ds), '_noavg.pdf', sep = ''),
         device = 'pdf',
         dpi = 600,
         width = 350,
         height = 800,
         units = 'mm',
         limitsize = F)
  
}

# use Desai's formula to quantify expected slopes from deltaF and eps_ij measurements
expected_slopes <- data.frame(dataset = character(0),
                              species = character(0),
                              expected_slope = numeric(0))
for (ds in unique(params$dataset)) {
  for (i in unique(params$species_i[params$dataset == ds])) {
    
    params_i <- params[params$dataset == ds & params$species_i == i, ]
    w_ij <- params_i$deltaF_j^2 / sum(params_i$deltaF_j^2)
    b_ij <- params_i$eps_ij / params_i$deltaF_j_Bi
    
    expected_slope <- sum(params_i$eps_ij * params_i$deltaF_j / sum(params_i$deltaF_j^2))
    
    if (!is.na(expected_slope)) {
      expected_slopes <- rbind(expected_slopes,
                               data.frame(dataset = ds,
                                          species = i,
                                          expected_slope = expected_slope))
    }
    
  }
}

# merge data frames and plot
expected_slopes$dataset[expected_slopes$dataset == basename(files)[6]] <- basename(files)[4]
fees$dataset[fees$dataset == basename(files)[6]] <- basename(files)[4]

plot_this <- merge(fees, expected_slopes, by = c('dataset', 'species'))

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

# plot slopes
lims <- c(min(c(plot_this$slope, plot_this$expected_slope)),
          max(c(plot_this$slope, plot_this$expected_slope)))

err <- sd(plot_this$slope - plot_this$expected_slope)

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
  scale_x_continuous(name = 'FEE slope\nexpected from theory',
                     limits = lims) +
  scale_y_continuous(name = 'Empirical FEE slope',
                     limits = lims) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank())

ggsave(filename = '../plots/slopes_predictions_noavg.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)


# write params table for comparison with phi_ij
write.table(params, file = 'params_phi.txt', quote = F, sep = '\t', row.names = F)

mycor <- do.call(rbind,
                 lapply(unique(plot_this$dataset),
                        FUN = function(ds) data.frame(dataset = ds,
                                                      spearman_cor = cor(plot_this$slope[plot_this$dataset == ds],
                                                                         plot_this$expected_slope[plot_this$dataset == ds],
                                                                         method = 'spearman'))))
write.table(mycor, file = 'corr_eps.txt', quote = F, sep = '\t', row.names = F)
