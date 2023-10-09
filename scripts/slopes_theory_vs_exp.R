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

# save FEE params
write.table(fees, file = './empirical_FEEs.txt', sep = '\t', row.names = F, quote = F)

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

plot_this$product <- abs(plot_this$deltaF_j_Bi * plot_this$eps_ij)
plot_this$product_sign <- as.character(sign(plot_this$deltaF_j_Bi * plot_this$eps_ij))

plot_this <- plot_this[!is.na(plot_this$product), ]

# plot starch dataset (for main fig)
ggplot(plot_this[plot_this$dataset == 'Bacterial starch hydrolysis', ],
       aes(x = deltaF_j_Bi, y = eps_ij, size = product, color = product_sign)) +
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
  scale_x_continuous(name = expression(paste(symbol("\341"), Delta, italic(F)[italic(j)], symbol("\361")[italic(B)(italic(i))], sep = '')),
                     breaks = pretty_breaks(n = 2),
                     expand = c(0.3, 0.3)) +
  scale_y_continuous(name = expression(paste(symbol("\341"), italic(epsilon)[italic(ij)], symbol("\361"), sep = '')),
                     breaks = pretty_breaks(n = 3),
                     expand = c(0.3, 0.3)) +
  scale_size(range = c(3, 8)) +
  scale_color_manual(values = c('firebrick1', 'deepskyblue')) +
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

ggsave(filename = '../plots/interaction_map_starch.pdf',
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
         aes(x = deltaF_j_Bi, y = eps_ij, size = product, color = product_sign)) +
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
    scale_x_continuous(name = expression(paste(symbol("\341"), Delta, italic(F)[italic(j)], symbol("\361")[italic(B)(italic(i))], sep = '')),
                       breaks = pretty_breaks(n = 2),
                       expand = c(0.3, 0.3)) +
    scale_y_continuous(name = expression(paste(symbol("\341"), italic(epsilon)[italic(ij)], symbol("\361"), sep = '')),
                       breaks = pretty_breaks(n = 3),
                       expand = c(0.3, 0.3)) +
    scale_size(range = c(3, 8)) +
    scale_color_manual(values = c('firebrick1', 'deepskyblue')) +
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
  
  ggsave(filename = paste('../plots/interaction_map_', gsub(' ', '', ds), '.pdf', sep = ''),
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
for (ds in names(data)) {
  
  ge_data_i <- makeGEdata(matrix2string(data[[ds]]))
  
  for (i in unique(ge_data_i$knock_in)) {
    
    avg_deltaF_i <- mean(ge_data_i$d_f[ge_data_i$knock_in == i])
    avg_F_Bi <- mean(ge_data_i$background_f[ge_data_i$knock_in == i])
    
    b_i <- expected_slopes$expected_slope[expected_slopes$dataset == ds & expected_slopes$species == i]
    
    expected_intercept <- avg_deltaF_i - b_i*avg_F_Bi
    
    if (length(expected_intercept)) {
      if (!is.na(expected_intercept)) {
        expected_intercepts <- rbind(expected_intercepts,
                                     data.frame(dataset = ds,
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
  } else {
    avg_deltaF_i <- mean(ge_data_i$d_f[ge_data_i$knock_in == 'P' & !grepl('T', ge_data_i$background)])
    avg_F_Bi <- mean(ge_data_i$background_f[ge_data_i$knock_in == 'P' & !grepl('T', ge_data_i$background)])
  }
  
  b_i <- expected_slopes$expected_slope[expected_slopes$dataset == unique(params$dataset)[ds] & expected_slopes$species == paste('P', thurin, sep = '_')]
  
  expected_intercept <- avg_deltaF_i - b_i*avg_F_Bi
  
  expected_intercepts <- rbind(expected_intercepts,
                               data.frame(dataset = unique(params$dataset)[ds],
                                          species = paste('P', thurin, sep = '_'),
                                          expected_intercept = expected_intercept))
  
}

# quantify effective interactions
sumF2 <- data.frame(dataset = character(0),
                    species_i = character(0),
                    sumF2 = numeric(0))

for (ds in unique(params$dataset)) {
  for (sp_i in unique(params$species_i[params$dataset == ds])) {
    
    allF <- params$deltaF_j_Bi[params$dataset == ds & params$species_i == sp_i]
    sumF2 <- rbind(sumF2,
                   data.frame(dataset = ds,
                              species_i = sp_i,
                              sumF2 = sum(allF^2)))
    
  }
}

params <- merge(params, sumF2, all = T)
params$eff_inter <- params$eps_ij * params$deltaF_j_Bi / params$sumF2

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

mycor <- do.call(rbind,
                 lapply(unique(plot_this$dataset),
                        FUN = function(ds) data.frame(dataset = ds,
                                                      spearman_cor = cor(plot_this$slope[plot_this$dataset == ds],
                                                                         plot_this$expected_slope[plot_this$dataset == ds],
                                                                         method = 'spearman'))))

# standardize intercepts so they can be plotted in a common scale
plot_this$intercept.std <- plot_this$intercept
plot_this$expected_intercept.std <- plot_this$expected_intercept
for (ds in unique(plot_this$dataset)) {
  
  avg <- mean(plot_this$intercept[plot_this$dataset == ds])
  std <- sd(plot_this$intercept[plot_this$dataset == ds])
  
  plot_this$intercept.std[plot_this$dataset == ds] <- (plot_this$intercept[plot_this$dataset == ds] - avg)/std
  plot_this$expected_intercept.std[plot_this$dataset == ds] <- (plot_this$expected_intercept[plot_this$dataset == ds] - avg)/std
  
}

# standardize slopes as well
plot_this$slope.std <- plot_this$slope
plot_this$expected_slope.std <- plot_this$expected_slope
for (ds in unique(plot_this$dataset)) {
  
  avg <- mean(plot_this$slope[plot_this$dataset == ds])
  std <- sd(plot_this$slope[plot_this$dataset == ds])
  
  plot_this$slope.std[plot_this$dataset == ds] <- (plot_this$slope[plot_this$dataset == ds] - avg)/std
  plot_this$expected_slope.std[plot_this$dataset == ds] <- (plot_this$expected_slope[plot_this$dataset == ds] - avg)/std
  
}

# plot slopes (standardized)
lims <- c(min(c(plot_this$slope.std, plot_this$expected_slope.std)),
          max(c(plot_this$slope.std, plot_this$expected_slope.std)))

err <- sd(plot_this$slope.std - plot_this$expected_slope.std)

ggplot(plot_this[!(plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P'), ],
       aes(x = expected_slope.std, y = slope.std, color = dataset, shape = type)) +
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
  scale_x_continuous(name = 'Standardized FEE slope\nexpected from pervasive pairwise epistasis',
                     limits = lims) +
  scale_y_continuous(name = 'Standardized observed FEE slope',
                     limits = lims) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank())

ggsave(filename = '../plots/slopes_predictions_std.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# plot slopes (not standardized)
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

ggsave(filename = '../plots/slopes_predictions.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# plot intercepts (standardized)
lims <- c(min(c(plot_this$intercept.std, plot_this$expected_intercept.std)),
          max(c(plot_this$intercept.std, plot_this$expected_intercept.std)))

err <- sd(plot_this$intercept.std - plot_this$expected_intercept.std)

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

ggsave(filename = '../plots/intercepts_predictions_std.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# plot intercepts (not standardized)
lims <- c(min(c(plot_this$intercept, plot_this$expected_intercept)),
          max(c(plot_this$intercept, plot_this$expected_intercept)))

err <- sd(plot_this$intercept - plot_this$expected_intercept)

ggplot(plot_this[!(plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species == 'P'), ],
       aes(x = expected_intercept, y = intercept, color = dataset, shape = type)) +
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
  scale_x_continuous(name = 'FEE intercept\nexpected from pervasive pairwise epistasis',
                     limits = lims) +
  scale_y_continuous(name = 'Observed FEE intercept',
                     limits = lims) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank())

ggsave(filename = '../plots/intercepts_predictions.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# for the starch dataset, plot FEEs with expected vs. observed fits
ge_data_i <- makeGEdata(matrix2string(data[[1]]))
sp <- unique(ge_data_i$knock_in)

# attach polymyxa branches
ge_data_i$branch <- T

attach_this <- ge_data_i[ge_data_i$knock_in == 'P', ]
attach_this$branch <- !grepl('T', attach_this$background)
attach_this$knock_in <- 'P_0'

ge_data_i <- rbind(ge_data_i, attach_this)

attach_this <- ge_data_i[ge_data_i$knock_in == 'P', ]
attach_this$branch <- grepl('T', attach_this$background)
attach_this$knock_in <- 'P_1'

ge_data_i <- rbind(ge_data_i, attach_this)

exp_i <- plot_this[plot_this$dataset == 'Bacterial starch hydrolysis', c('species', 'slope', 'intercept', 'expected_slope', 'expected_intercept')]
colnames(exp_i)[1] <- 'knock_in'

ge_data_i$knock_in <- sp_names[ge_data_i$knock_in]
exp_i$knock_in <- sp_names[exp_i$knock_in]

ggplot(ge_data_i, aes(x = background_f, y = d_f, group = knock_in, alpha = branch)) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3,
             shape = 16) +
  geom_abline(data = exp_i,
              aes(slope = slope, intercept = intercept),
              linetype = 'solid',
              color = 'black',
              size = 1) +
  geom_abline(data = exp_i,
              aes(slope = expected_slope, intercept = expected_intercept),
              linetype = 'dashed',
              color = '#d6d62d',
              size = 1) +
  scale_x_continuous(breaks = pretty_breaks(n = 2),
                     name = 'Function of ecological background [a.u.]') +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = expression(paste(Delta, italic(F), ' [a.u.]', sep = ''))) +
  scale_alpha_manual(values = c(0.1, 1)) +
  facet_wrap(~knock_in,
             ncol = 6) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 12,
                                  vjust = 0),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = 'none')

ggsave(filename = '../plots/starch_predicted_fees.pdf',
       device = 'pdf',
       dpi = 600,
       width = 350,
       height = 800,
       units = 'mm',
       limitsize = F)

# plot three examples:
plot_this <- plot_this[plot_this$species %in% c('fa.nat', 'SL106', 'P_1'), ]
plot_this$species <- sp_names[plot_this$species]
plot_this$species[1] <- 'P. polymyxa'

ge_data_i <- makeGEdata(matrix2string(data[[1]]))
ge_data_i <- ge_data_i[ge_data_i$knock_in == 'P', ]
ge_data_i$branch <- grepl('T', ge_data_i$background)

ge_data_i <- rbind(ge_data_i, cbind(ge_data[[4]], branch = F))
ge_data_i <- rbind(ge_data_i, cbind(ge_data[[5]], branch = F))
ge_data_i <- ge_data_i[ge_data_i$knock_in %in% c('fa.nat', 'SL106', 'P'), ]
ge_data_i$knock_in <- sp_names[ge_data_i$knock_in]

colnames(ge_data_i)[2] <- 'species'
ge_data_i$species <- factor(ge_data_i$species,
                            levels = c('P. polymyxa',
                                       'Sphingoterrabacterium sp.',
                                       'L. capitata'))

lims <- data.frame(species = rep(levels(ge_data_i$species), 2),
                   background_f = c(20, 1.5, 5.5, 0, 0.5, 0),
                   d_f = c(32, 1.2, 12, -8, -0.25, -3))
lims$species <- factor(lims$species, levels = levels(ge_data_i$species))
lims <- rbind(cbind(lims, branch = T),
              cbind(lims, branch = F))

plot_this$species <- factor(plot_this$species, levels = levels(ge_data_i$species))

ggplot(ge_data_i, aes(x = background_f, y = d_f, group = branch)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             shape = 1,
             cex = 3) +
  geom_point(data = lims,
             color = 'white',
             shape = 3) +
  geom_smooth(formula = y ~ x,
              method = 'lm',
              color = 'black',
              fullrange = F,
              se = F) +
  # geom_abline(data = plot_this,
  #             aes(slope = slope, intercept = intercept)) +
  geom_abline(data = plot_this,
              aes(slope = expected_slope, intercept = expected_intercept),
              linetype = 'dashed') +
  facet_wrap(~species,
             ncol = 1,
             scales = 'free') +
  scale_x_continuous(name = expression(paste(italic(F), ' (background) [a.u.]', sep = '')),
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = expression(paste(Delta, italic(F), ' [a.u.]', sep = '')),
                     breaks = pretty_breaks(n = 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 16,
                                  hjust = 0),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/fees_expvstheory_samples.pdf',
       device = 'pdf',
       dpi = 600,
       width = 100,
       height = 150,
       units = 'mm',
       limitsize = F)

# plot effective interactions
plot_this <- params[params$species_i %in% c('fa.nat', 'SL106', 'P_1'), ]
plot_this$species_i <- sp_names[plot_this$species_i]
plot_this$species_j <- sp_names[plot_this$species_j]
plot_this$species_i[plot_this$species_i == 'P. polymyxa\n(B. thuringiensis\nin background)'] <- 'P. polymyxa'

plot_this$rank <- NA
for (sp_i in unique(plot_this$species_i)) {
  
  plot_this$rank[plot_this$species_i == sp_i] <- 6 - rank(plot_this$eff_inter[plot_this$species_i == sp_i])
  
}

plot_this$species_i <- factor(plot_this$species_i, levels = levels(ge_data_i$species))
plot_this$ytxt <- pmax(plot_this$eff_inter, 0)

ggplot(plot_this,
       aes(x = rank, y = eff_inter, fill = eff_inter > 0, label = species_j)) +
  geom_bar(stat = 'identity',
           width = 0.5) +
  geom_hline(yintercept = 0) +
  geom_text(aes(y = ytxt),
            hjust = -0.1,
            fontface = 'italic',
            size = 5) +
  facet_wrap(~species_i,
             ncol = 1,
             scales = 'free') +
  scale_x_continuous(limits = c(0.5, 5.5)) +
  scale_y_continuous(name = expression(paste('Effective interaction, ', tilde(italic(epsilon))[italic(ij)], sep = '')),
                     breaks = pretty_breaks(n = 2)) +
  scale_fill_manual(values = c('firebrick1', 'deepskyblue')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 16,
                                  hjust = 0),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  coord_flip(clip = 'off') +
  annotate("segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=0.5)

ggsave(filename = '../plots/interMap_samples.pdf',
       device = 'pdf',
       dpi = 600,
       width = 100,
       height = 150,
       units = 'mm',
       limitsize = F)

# plot polymyxa's effective interactions
plot_this <- params[params$species_i %in% c('P', 'P_0', 'P_1'), ]
plot_this$species_i <- sp_names[plot_this$species_i]
plot_this$species_j <- sp_names[plot_this$species_j]

plot_this$rank <- NA
for (sp_i in unique(plot_this$species_i)) {
  
  plot_this$rank[plot_this$species_i == sp_i] <- 6 - rank(plot_this$eff_inter[plot_this$species_i == sp_i])
  
}

plot_this$species_i <- factor(plot_this$species_i, levels = unique(plot_this$species_i))
plot_this$ytxt <- pmax(plot_this$eff_inter, 0)

lims <- data.frame(species_i = rep(unique(plot_this$species_i), 2),
                   eff_inter = c(-0.67, -6.7, -6.7, 0.67, 6.7, 6.7),
                   rank = rep(2.5, 6),
                   species_j = 'scale_limits')

ggplot(plot_this,
       aes(x = rank, y = eff_inter, fill = eff_inter > 0, label = species_j)) +
  geom_bar(stat = 'identity',
           width = 0.5) +
  geom_hline(yintercept = 0) +
  geom_text(aes(y = ytxt),
            hjust = -0.1,
            fontface = 'italic',
            size = 5) +
  facet_wrap(~species_i,
             ncol = 1,
             scales = 'free',
             labeller = labeller(species_i = setNames(c('P. polymyxa', 'P. polymyxa (no T)', 'P. polymyxa (T)'),
                                                      unique(plot_this$species_i)))) +
  geom_blank(data = lims) +
  scale_x_continuous(limits = c(0.5, 5.5)) +
  scale_y_continuous(name = expression(paste('Effective interaction, ', tilde(italic(epsilon))[italic(ij)], sep = '')),
                     breaks = pretty_breaks(n = 3)) +
  scale_fill_manual(values = c('firebrick1', 'deepskyblue')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 16,
                                  hjust = 0),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  coord_flip(clip = 'off') +
  annotate("segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=0.5)

ggsave(filename = '../plots/interMap_polymyxa.pdf',
       device = 'pdf',
       dpi = 600,
       width = 100,
       height = 150,
       units = 'mm',
       limitsize = F)

######################################################
### functional landscape
######################################################

# this chunk of code is borrowed from:
# Juan Diaz-Colunga, Abigail Skwara, Karna Gowda, Ramon Diaz-Uriarte,
# Mikhail Tikhonov, Djordje Bajic & Alvaro Sanchez (2022). Global epistasis
# on fitness landscapes. arXiv. https://arxiv.org/abs/2210.03677

# tunable parameters
n_mut <- 6
mycolors <- c('#939598', '#be1e2d', '#85c441', '#d68f28', '#415ba9', '#a96cad')
mycolors <- c('#939598', '#d68f28', '#415ba9', '#a96cad')

# make genotype names
genots <- lapply(0:n_mut, FUN = function(i) t(combn(n_mut, i)))
genots <- lapply(genots, FUN = function(x) sapply(1:nrow(x),
                                                  FUN = function(i) paste(x[i, ], collapse = ',')))
genots <- unlist(genots)

# make edges of fitness graph
nMut <- function(genot) sapply(genot,
                               FUN = function(genot_i) length(strsplit(genot_i, split = ',')[[1]]))

isDescendant <- function(this_genot, of_this_genot) {
  
  this_genot <- strsplit(this_genot, split = ',')[[1]]
  of_this_genot <- strsplit(of_this_genot, split = ',')[[1]]
  
  return(all(of_this_genot %in% this_genot))
  
}

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

# plot landscape (wrapper function)
plotGraph <- function(landscape, save.plot = F) {
  
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
  
  myplot <-
    ggplot(df, aes(x = nmut, y = f, group = edge_id, color = color)) +
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

landscape <- data.frame(genot = genots)
data_i <- matrix2string(data[['amyl_Sanchez-Gorostiaga2019.csv']])

sp_names <- unique(params$species_i[params$dataset == 'amyl_Sanchez-Gorostiaga2019.csv'])
sp_names <- sp_names[!grepl('P_', sp_names)]

for (i in 1:6) data_i$community <- gsub(sp_names[i], i, data_i$community)
landscape <- merge(landscape, data_i, by.x = 'genot', by.y = 'community', all = T)
colnames(landscape)[2] <- 'f'

plotGraph(landscape) # this spits out a warning because the landscape is combinatorially incomplete, but the plot is produced

ggsave(file = paste('../plots/fitnessGraph-full_.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')

# inset
df <- matrix2string(data[['amyl_Sanchez-Gorostiaga2019.csv']])
colnames(df) <- c('genot', 'f')
df <- df[df$genot %in% c('C,M', 'C,M,P', 'C,M,T', 'C,M,P,T'), ]
df$nmut <- sapply(df$genot, FUN = function(x) length(strsplit(x, split = ',')[[1]]))

ggplot(df, aes(x = nmut, y = f)) +
  geom_point(cex = 3) +
  scale_y_continuous(name = 'Function [a.u.]',
                     limits = c(min(df$f), max(c(df$f[df$nmut == max(df$nmut)],
                                                 sum(df$f[df$nmut == (max(df$nmut) - 1)]) - df$f[df$nmut == min(df$nmut)]))),
                     expand = rep(0.05, 2),
                     breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(expand = rep(0.1, 2)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, vjust = 0),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.5)

ggsave(file = paste('../plots/fitnessGraph_inset.pdf', sep = ''),
       device = 'pdf',
       width = 100,
       height = 100,
       units = 'mm')




# write params table for comparison with phi_ij
write.table(params, file = 'params_eps.txt', quote = F, sep = '\t', row.names = F)

