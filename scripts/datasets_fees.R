# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(cowplot)

# load data sets
files <- list.files('../data_sets', full.names = T)
files <- c(files, '../pyoverdine_data/training_set.csv')
files <- c(files, '../genetic_data_sets/Khan_fitness.csv')
data <- lapply(files, FUN = function(file) read.csv(file))

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

# we first make dF vs F_background plots for our data (pyoverdine secretion)
pyo_data <- data[[6]]
pyo_data <- lapply(1:3,
                   FUN = function(i) matrix2string(pyo_data[1:(nrow(pyo_data) - 1), c(1:8, 8+i)]))

ge_data <- lapply(pyo_data, FUN = makeGEdata)
ge_data <- merge(merge(ge_data[[1]], ge_data[[2]],
                       by = c('background', 'knock_in'),
                       suffixes = c('.rep1', '.rep2'),
                       all = T),
                 ge_data[[3]],
                 by = c('background', 'knock_in'),
                 suffixes = c('', '.rep3'),
                 all = T)
colnames(ge_data)[7:8] <- c('background_f.rep3', 'd_f.rep3')

ge_data$background_f.mean <- as.numeric(sapply(1:nrow(ge_data),
                                               FUN = function(i) mean(as.numeric(ge_data[i, paste('background_f.rep', 1:3, sep = '')]))))
ge_data$background_f.sd <- as.numeric(sapply(1:nrow(ge_data),
                                             FUN = function(i) sd(as.numeric(ge_data[i, paste('background_f.rep', 1:3, sep = '')]))))
ge_data$d_f.mean <- as.numeric(sapply(1:nrow(ge_data),
                                      FUN = function(i) mean(as.numeric(ge_data[i, paste('d_f.rep', 1:3, sep = '')]))))
ge_data$d_f.sd <- as.numeric(sapply(1:nrow(ge_data),
                                    FUN = function(i) sd(as.numeric(ge_data[i, paste('d_f.rep', 1:3, sep = '')]))))

ge_data$knock_in <- sp_names[[6]][ge_data$knock_in]
ge_data$knock_in <- factor(ge_data$knock_in,
                           levels = c("Pseudomonas sp. 01",
                                      "Pseudomonas sp. 02",
                                      "Pseudomonas sp. 03",
                                      "Pseudomonas sp. 04",
                                      "Pseudomonas sp. 05",
                                      "Enterobacter sp.",
                                      "Raoultella sp.",
                                      "Klebsiella sp.")) # ordered from higher to lower function in monoculture

# plot FEEs
ggplot(ge_data, aes(x = background_f.mean, xmin = background_f.mean - background_f.sd, xmax = background_f.mean + background_f.sd,
                    y = d_f.mean, ymin = d_f.mean - d_f.sd, ymax = d_f.mean + d_f.sd,
                    group = knock_in)) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 2,
             shape = 16) +
  geom_errorbar(alpha = 0.25) +
  geom_errorbarh(alpha = 0.25) +
  geom_smooth(method = 'lm',
              formula = y~x,
              color = 'firebrick1',
              se = F,
              fullrange = T) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = 'Function of ecological background [uM]') +
  scale_y_continuous(breaks = pretty_breaks(n = 3),
                     name = expression(paste(Delta, italic(F), ' [uM]', sep = ''))) +
  facet_wrap(~knock_in,
             nrow = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 16),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/fees_pyoverdine.pdf',
       device = 'pdf',
       dpi = 600,
       width = 160,
       height = 100,
       units = 'mm',
       limitsize = F)

# also plot global epistasis patterns of Khan et al. dataset (for reference, this goes in fig. 1)
ge_data <- makeGEdata(matrix2string(data[[7]]))
ge_data$knock_in <- sp_names[[7]][ge_data$knock_in]

ge_data <- ge_data[ge_data$knock_in != 'topA', ] # we only plot 4 mutations
ge_data$knock_in <- factor(ge_data$knock_in,
                           levels = c('glmUS', 'pykF', 'rbs', 'spoT'))

ggplot(ge_data, aes(x = background_f, y = d_f, group = knock_in)) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 3,
             shape = 16) +
  geom_smooth(method = 'lm',
              formula = y~x,
              color = 'firebrick1',
              se = F,
              fullrange = T) +
  scale_x_continuous(breaks = pretty_breaks(n = 2),
                     name = expression(paste('Relative fitness of genetic background, ', italic(f)[B], sep = ''))) +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = expression(paste(Delta, ' fitness', sep = ''))) +
  facet_wrap(~knock_in,
             nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 16),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/fees_khan.pdf',
       device = 'pdf',
       dpi = 600,
       width = 160,
       height = 100,
       units = 'mm',
       limitsize = F)

# make plots of dF vs. F_background (all other datasets, this goes in fig. 2)
g <- vector(mode = 'list', length = length(data))

for (i in 1:5) {
  
  # if there are multiple measurements of a same community, average them
  colnames(data[[i]])[ncol(data[[i]])] <- 'fun'
  data[[i]] <- aggregate(formula = fun ~ .,
                         data = data[[i]],
                         FUN = mean)
  
  # for the phytoplankton biomass dataset (Ghedini et al., scale functions by 1e-4 for easier readability)
  if (i == 3) data[[i]][, ncol(data[[i]])] <- data[[i]][, ncol(data[[i]])]/1e4
  
  # extract F vs. dF data from combinatorial assemblages
  ge_data <- makeGEdata(matrix2string(data[[i]]))
  
  # remove species for which there are not at least 5 data points (most data sets are combinatorially incomplete)
  npoints <- sapply(unique(ge_data$knock_in),
                    FUN = function(sp) sum(ge_data$knock_in == sp))
  valid_sp <- names(npoints)[npoints >= 5]
  ge_data <- ge_data[ge_data$knock_in %in% valid_sp, ]
  
  # group for linear fits (group by species unless branching is observed, i.e. P. polymyxa when B. thruingiensis is present in amyl data)
  ge_data$group <- ge_data$knock_in
  
  if (any(ge_data$knock_in == 'P')) ge_data$group[ge_data$knock_in == 'P' & grepl('T', ge_data$background)] <- 'P.T'
  
  # make linear fits and get slopes
  slopes <- sapply(unique(ge_data$group),
                   FUN = function(g) as.numeric(lm(ge_data$d_f[ge_data$group == g] ~
                                                     ge_data$background_f[ge_data$group == g])$coefficients[2]))
  slopes[slopes > 0.8] <- 0.8 # 'shrink' slope values from the tails of the distribution to avoid NAs in color scale
  slopes[slopes < -0.8] <- -0.8
  ge_data$slope <- slopes[ge_data$group]
  
  # manually set axis limits for clear visualization
  dx <- (max(ge_data$background_f) - min(ge_data$background_f))/20
  dy <- (max(ge_data$d_f) - min(ge_data$d_f))/20
  
  dx <- c(min(ge_data$background_f) - dx, max(ge_data$background_f) + dx)
  dy <- c(min(ge_data$d_f) - dy, max(ge_data$d_f) + dy)
  
  # full species names
  ge_data$knock_in <- sp_names[[i]][ge_data$knock_in]
  
  # if there are less than 25 species (maximum across our data sets), add 'ghost' panels to mantain axis proportions
  if (length(unique(ge_data$knock_in)) < 25) {
    
    add_panels <- 25 - length(unique(ge_data$knock_in))
    
    ge_data <- rbind(ge_data,
                     data.frame(background = NA,
                                knock_in = paste('remove_this_panel', 1:add_panels, sep = '.'),
                                background_f = NA,
                                d_f = NA,
                                group = NA,
                                slope = NA))
    
  }
  ge_data$knock_in <- factor(ge_data$knock_in, levels = unique(ge_data$knock_in))
  
  # plot
  g[[i]] <- 
    ggplot(ge_data, aes(x = background_f, y = d_f, group = group, color = slope)) +
    geom_abline(slope = 0,
                intercept = 0,
                linetype = 'dashed') +
    geom_point(color = 'black',
               shape = 20) +
    geom_smooth(method = 'lm',
                formula = y~x,
                se = FALSE,
                fullrange = TRUE) +
    scale_x_continuous(name = '',
                       limits = dx,
                       breaks = pretty_breaks(n = 2)) +
    scale_y_continuous(name = 'dF',
                       limits = dy,
                       breaks = pretty_breaks(n = 2)) +
    facet_wrap(~ knock_in,
               ncol = 7) +
    scale_color_gradientn(colors = c('firebrick1',
                                     'firebrick1',
                                     'black',
                                     'deepskyblue',
                                     'deepskyblue'),
                          breaks = c(-0.8, -0.5, 0, 0.5, 0.8),
                          limits = c(-0.8, 0.8),
                          na.value = 'deepskyblue') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
  
}

myplot <- grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 5)

ggsave(myplot,
       filename = '../plots/datasets_fees.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 500,
       units = 'mm',
       limitsize = F)

# quality of fits
# for species i, quality is quantified as 1 - sigma_fee^2/sigma^2, where:
#   sigma^2 is the variance of the distribution of functional effects across all species in the dataset
#   sigma_fee^2 is the variance of the distribution of residuals for the linear fit of species i

sigma <- data.frame(dataset = character(0),
                    species = character(0),
                    sigma_fee = numeric(0),
                    sigma = numeric(0))

for (i in 1:6) {
  
  ge_data <- makeGEdata(matrix2string(data[[i]]))
  
  ge_data$group <- ge_data$knock_in
  if (any(ge_data$knock_in == 'P')) ge_data$group[ge_data$knock_in == 'P' & grepl('T', ge_data$background)] <- 'P.T'
  
  # make linear fits
  lfits <- lapply(unique(ge_data$group),
                  FUN = function(grp) lm(formula = d_f ~ background_f,
                                         data = ge_data[ge_data$group == grp, ])$coefficients)
  names(lfits) <- unique(ge_data$group)
  
  for (grp in unique(ge_data$group)) {
    yfit <- lfits[[grp]][1] + lfits[[grp]][2] * ge_data$background_f[ge_data$group == grp]
    sigma <- rbind(sigma,
                   data.frame(dataset = basename(files[i]),
                              species = grp,
                              sigma_fee = sd(yfit - ge_data$d_f[ge_data$group == grp]),
                              sigma = sd(ge_data$d_f)))
  }
  
}

sigma$quality <- 1 - sigma$sigma_fee^2/sigma$sigma^2
sigma$dataset <- factor(sigma$dataset,
                        levels = c("training_set.csv",
                                   "plant-biommass_Kuebbing2016_all.csv",
                                   "phytoplankton-biomass_Ghedini2022.csv",
                                   "xylose_Langenheder2010.csv",
                                   "amyl_Sanchez-Gorostiaga2019.csv",
                                   "butyrate_Clark2021.csv"))

# plot
myplot1 <- 
  ggplot(sigma, aes(x = 0, y = quality)) +
  geom_jitter(size = 3) +
  facet_wrap(~dataset, nrow = 1) +
  scale_x_continuous(name = '',
                     limits = c(-0.75, 0.75)) +
  scale_y_continuous(name = expression(1 - italic(sigma)[fee]^2 / sigma^2),
                     limits = c(0, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 3,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.5) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5)

myplot2 <-
  ggplot(sigma, aes(x = quality)) +
  geom_histogram(bins = 20) +
  coord_flip() +
  scale_x_continuous(name = '') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 3,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=0.5)

r <- 0.81
myplot <- plot_grid(myplot1, myplot2, rel_widths = c(r, 1-r))

print(myplot)
ggsave(myplot,
       filename = '../plots/datasets_fees_quality.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)















# # distribution of R squared
# ge_data <- lapply(data, FUN = function(x) makeGEdata(matrix2string(x)))
# data_sets <- gsub('.csv', '', basename(files))
# ge_data <- lapply(1:length(ge_data), FUN = function(i) cbind(data_set = data_sets[i], ge_data[[i]]))
# ge_data <- do.call(rbind, ge_data)
# ge_data$knock_in[ge_data$knock_in == 'P' & containsSpecies('T', ge_data$background)] <- 'P-1'
# 
# rsq <- data.frame(data_set = character(0),
#                   species = character(0),
#                   rsq = numeric(0),
#                   slope = numeric(0))
# for (ds in data_sets) {
#   for (sp in unique(ge_data$knock_in[ge_data$data_set == ds])) {
#     
#     n <- ge_data$data_set == ds & ge_data$knock_in == sp
#     rsq <- rbind(rsq,
#                  data.frame(data_set = ds,
#                             species = sp,
#                             rsq = cor(ge_data$background_f[n], ge_data$d_f[n])^2,
#                             slope = lm(d_f ~ background_f, data = ge_data[n, ])$coefficients[2]))
#     
#   }
# }
# rsq$species[rsq$species == 'P-1'] <- 'P'
# rsq$data_set <- factor(rsq$data_set, levels = data_sets[c(4, 3, 5, 1, 2)])
# 
# myplot1 <- 
#   ggplot(rsq, aes(x = 0, y = rsq)) +
#   geom_jitter(size = 3) +
#   facet_wrap(~data_set, nrow = 1) +
#   scale_x_continuous(name = '',
#                      limits = c(-0.75, 0.75)) +
#   scale_y_continuous(name = expression(paste(italic(R)^2, ' of FEE')),
#                      limits = c(0, 1)) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(face = 'italic',
#                                   size = 10),
#         aspect.ratio = 3,
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 18),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         legend.position = 'none') +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
#   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
#   annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, size=0.5) +
#   annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5)
# 
# myplot2 <-
#   ggplot(rsq, aes(x = rsq)) +
#   geom_histogram(bins = 20) +
#   coord_flip() +
#   scale_x_continuous(name = expression(paste(italic(R)^2, ' of FEE'))) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(face = 'italic',
#                                   size = 10),
#         aspect.ratio = 3,
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 18),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         legend.position = 'none') +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=0.5)
# 
# r <- 0.75
# myplot <- plot_grid(myplot1, myplot2, rel_widths = c(r, 1-r))
# 
# print(myplot)
# ggsave(myplot,
#        filename = '../plots/figS0.5.pdf',
#        device = 'pdf',
#        dpi = 600,
#        width = 200,
#        height = 100,
#        units = 'mm',
#        limitsize = F)
