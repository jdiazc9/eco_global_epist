rm(list = ls())
library(tidyr)
library(ggplot2)

files <- list.files('../data_sets', full.names = T)

# load params files
phi <- read.table('params_phi.txt', header = T, sep = '\t')
eps <- read.table('params_eps.txt', header = T, sep = '\t')

# calculate phi (effective interactions estimated just from monocultures and pairs data)
phi_sumF2 <- aggregate(deltaF_j ~ dataset + species_i,
                       data = phi,
                       FUN = function (x) sum(x^2))
phi <- merge(phi, phi_sumF2, by = c('dataset', 'species_i'),
             suffixes = c('', '.sum2'))
phi$phi <- phi$eps_ij * phi$deltaF_j / phi$deltaF_j.sum2

phi <- phi[, c('dataset', 'species_i', 'species_j', 'phi')]
eps <- eps[, c('dataset', 'species_i', 'species_j', 'eff_inter')]

eps$dataset <- gsub('_inv', '', eps$dataset)
phi$dataset <- gsub('_inv', '', phi$dataset)

# get estimated FEE slopes
est_slopes_eps <- aggregate(eff_inter ~ dataset + species_i,
                            data = eps,
                            FUN = function(x) sum(x, na.rm = T))
est_slopes_phi <- aggregate(phi ~ dataset + species_i,
                            data = phi,
                            FUN = function(x) sum(x, na.rm = T))

# get empirical FEE slopes and merge all
fees <- read.table('empirical_FEEs.txt', sep = '\t', header = T)

slopes <- merge(fees[, c('dataset', 'species', 'slope')],
                est_slopes_eps,
                by.x = c('dataset', 'species'),
                by.y = c('dataset', 'species_i'),
                all = T)
slopes <- merge(slopes, est_slopes_phi,
                by.x = c('dataset', 'species'),
                by.y = c('dataset', 'species_i'),
                all = T)
colnames(slopes) <- c('dataset', 'species', 'empirical_slope', 'slope_eps', 'slope_phi')

# clean data frame for plotting
plot_this <- merge(phi, eps, all = T)
#plot_this <- plot_this[!is.na(plot_this$phi) & !is.na(plot_this$eff_inter), ]

plot_this$dataset <- gsub('_inv', '', plot_this$dataset)
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
plot_this$species_j <- sp_names[plot_this$species_j]

# plot example: effective interactions for a particular species (phi vs. eps)
focal_sp <- 'B. cereus'
focal_sp <- 'P. pratense'
df <- plot_this[plot_this$species_i == focal_sp & !is.na(plot_this$species_i), ]
if(focal_sp %in% c('B. subtilis', 'B. cereus', 'B. thuringiensis', 'B. mojavensis', 'P. polymyxa')) df <- df[df$dataset == 'Bacterial starch hydrolysis', ]
df <- gather(df, param, value, phi:eff_inter, factor_key = TRUE)
df$param <- factor(df$param, levels = c('eff_inter', 'phi'))

df_sum <- aggregate(value ~ dataset + species_i + param + type,
                    data = df,
                    FUN = sum)
df_sum$species_j <- 'Sum of all'
df_sum <- df_sum[, c('dataset', 'species_i', 'species_j', 'type', 'param', 'value')]
df <- rbind(df, df_sum)
df$species_j <- factor(df$species_j,
                       levels = c(unique(df$species_j[df$species_j != 'Sum of all']), 'Sum of all'))

df$xpos <- setNames(1:length(unique(df$species_j)),
                    levels(df$species_j))[df$species_j]
Dx <- 0.5
df$xpos[df$xpos == max(df$xpos)] <- max(df$xpos) + Dx

ggplot(df, aes(x = xpos, y = value, fill = param)) +
  geom_bar(stat = 'identity',
           position = position_dodge(),
           width = 0.6) +
  geom_hline(yintercept = 0,
             color = 'black') + 
  geom_vline(xintercept = max(df$xpos) - 3*Dx/2,
             color = 'black',
             linetype = 'dashed') +
  scale_x_reverse(name = expression(Species~italic(j)),
                     breaks = unique(df$xpos),
                     labels = unique(df$species_j),
                  position = 'top') +
  scale_fill_manual(values = c("#d6d52b", "#473c85"),
                    labels = c(expression(tilde(italic(epsilon))[italic(ij)]),
                              expression(italic(phi)[italic(ij)]))) +
  coord_flip() +
  ggtitle(focal_sp) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 16),
        axis.text.y = element_text(size = 14,
                                   angle = 0,
                                   hjust = 1,
                                   face = 'italic'),
        axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = 'bottom',
        plot.title = element_text(size = 16,
                                  face = 'italic'),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = paste('../plots/compare_phi_eps_', focal_sp, '.pdf', sep = ''),
       device = 'pdf',
       dpi = 600,
       width = 125,
       height = 80,
       units = 'mm',
       limitsize = F)

# plot all species and datasets
mycolors <- setNames(c('black',
                       '#d6d62d',
                       '#66b666',
                       '#cb96c3',
                       '#d72027',
                       '#519ed7',
                       'black'),
                     levels(plot_this$dataset))

lims <- range(c(plot_this$phi, plot_this$eff_inter), na.rm = T)
ggplot(plot_this[!(plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species_i == 'P. polymyxa'), ],
       aes(x = phi, y = eff_inter, color = dataset, shape = type)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_point(cex = 3) +
  scale_color_manual(values = mycolors,
                     breaks = names(mycolors),
                     name = 'Data set') +
  scale_shape_manual(values = c(16, 3),
                     name = 'Type of data set') +
  scale_x_continuous(name = expression(italic(phi)[italic(ij)]),
                     limits = lims) +
  scale_y_continuous(name = expression(tilde(italic(epsilon))[italic(ij)]),
                     limits = lims) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank())

ggsave(filename = '../plots/compare_phi_eps_all.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# break down by dataset
mycor <- do.call(rbind,
                 lapply(unique(plot_this$dataset),
                        FUN = function(ds) {
                          
                          x <- plot_this$phi[plot_this$dataset == ds]
                          y <- plot_this$eff_inter[plot_this$dataset == ds]
                          
                          n <- !is.na(x) & !is.na(y)
                          x <- x[n]
                          y <- y[n]
                          
                          res <- y - x
                          R2_identity <- max(1 - sum(res^2)/sum((y - mean(y))^2), 0) # R squared of the y = x model
                          
                          return(data.frame(dataset = ds,
                                            cor = R2_identity))
                          
                        }))
mycor$cor <- sapply(mycor$cor,
                    FUN = function(x) paste('R2 = ', substr(as.character(x), 1, 4), sep = ''))

ranges <- do.call(rbind,
                  lapply(unique(plot_this$dataset),
                         FUN = function(ds){
                           
                           plot_this_i <- plot_this[plot_this$dataset == ds, ]
                           return(data.frame(dataset = ds,
                                             min = min(c(plot_this_i$phi, plot_this_i$eff_inter), na.rm = T),
                                             max = max(c(plot_this_i$phi, plot_this_i$eff_inter), na.rm = T)))
                           
                         }))
mycor <- merge(mycor, ranges, all = T)
mycor <- merge(mycor, unique(plot_this[, c('dataset', 'type')]), all = T)

ggplot(plot_this[!(plot_this$dataset == 'Bacterial starch hydrolysis' & plot_this$species_i == 'P. polymyxa'), ],
       aes(x = phi, y = eff_inter, color = dataset, shape = type)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_point(cex = 3) +
  geom_blank(aes(x = eff_inter, y = phi)) + # same x and y scales
  geom_text(data = mycor,
            aes(x = max, y = min, label = cor),
            color = 'black',
            hjust = 1,
            vjust = 0) +
  scale_color_manual(values = mycolors,
                     breaks = names(mycolors),
                     name = 'Data set') +
  scale_shape_manual(values = c(16, 3),
                     name = 'Type of data set') +
  scale_x_continuous(name = expression(italic(phi)[italic(ij)]),
                     breaks = pretty_breaks(n = 2)) +
  scale_y_continuous(name = expression(tilde(italic(epsilon))[italic(ij)]),
                     breaks = pretty_breaks(n = 2)) +
  facet_wrap(~ dataset, nrow = 1, scales = 'free') +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = 'none')

ggsave(filename = '../plots/compare_phi_eps_bydataset.pdf',
       device = 'pdf',
       dpi = 600,
       width = 300,
       height = 100,
       units = 'mm',
       limitsize = F)

# compare slopes estimated by eps and phi
slopes$dataset <- setNames(c('Bacterial starch hydrolysis',
                                'Bacterial butyrate secretion',
                                'Phytoplankton biomass',
                                'Above-ground plant biomass',
                                'Bacterial xylose oxidation',
                                'Bacterial pyoverdine secretion',
                                'E. coli fitness'),
                              c(basename(files)[1:5], 'pyo', 'khan'))[slopes$dataset]
slopes$dataset <- factor(slopes$dataset, levels = c('E. coli fitness',
                                                          'Above-ground plant biomass',
                                                          'Phytoplankton biomass',
                                                          'Bacterial xylose oxidation',
                                                          'Bacterial starch hydrolysis',
                                                          'Bacterial butyrate secretion',
                                                          'Bacterial pyoverdine secretion'))

slopes$type <- c('Ecological data set', 'Genetic data set')[1 + slopes$dataset %in% c('E. coli fitness')]
slopes$type <- factor(slopes$type, levels = c('Genetic data set',
                                                    'Ecological data set'))

lims <- range(c(slopes$empirical_slope, slopes$slope_eps))
ggplot(slopes[!(slopes$dataset == 'Bacterial starch hydrolysis' & slopes$species == 'P'), ],
       aes(x = slope_eps, y = empirical_slope, color = dataset, shape = type)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_point(cex = 3) +
  scale_color_manual(values = mycolors,
                     breaks = names(mycolors),
                     name = 'Data set') +
  scale_shape_manual(values = c(16, 3),
                     name = 'Type of data set') +
  scale_x_continuous(name = 'FEE slope\nexpected from averaged interactions',
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

ggsave(filename = '../plots/slopes_predictions_eps.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# break down by dataset
mycor <- do.call(rbind,
                 lapply(unique(slopes$dataset),
                               FUN = function(ds){
                                 
                                 slopes_i <- slopes[slopes$dataset == ds, ]
                                 
                                 x <- slopes_i$empirical_slope
                                 y_eps <- slopes_i$slope_eps
                                 y_phi <- slopes_i$slope_phi
                                 
                                 R2_eps <- max(1 - sum((y_eps - x)^2)/sum((y_eps - mean(y_eps))^2), 0)
                                 R2_phi <- max(1 - sum((y_phi - x)^2)/sum((y_phi - mean(y_phi))^2), 0)
                                 
                                 return(data.frame(dataset = ds,
                                                   cor_eps = R2_eps,
                                                   cor_phi = R2_phi))
                                 
                               }))
mycor$cor_eps <- sapply(mycor$cor_eps,
                        FUN = function(x) paste('R2 = ', substr(as.character(x), 1, 4), sep = ''))
mycor$cor_phi <- sapply(mycor$cor_phi,
                        FUN = function(x) paste('R2 = ', substr(as.character(x), 1, 4), sep = ''))

ranges <- do.call(rbind,
                  lapply(unique(slopes$dataset),
                         FUN = function(ds){
                           
                           slopes_i <- slopes[slopes$dataset == ds, ]
                           return(data.frame(dataset = ds,
                                             min_eps = min(c(slopes_i$empirical_slope, slopes_i$slope_eps)),
                                             max_eps = max(c(slopes_i$empirical_slope, slopes_i$slope_eps)),
                                             min_phi = min(c(slopes_i$empirical_slope, slopes_i$slope_phi)),
                                             max_phi = max(c(slopes_i$empirical_slope, slopes_i$slope_phi))))
                           
                         }))
mycor <- merge(mycor, ranges, all = T)
mycor <- merge(mycor, unique(slopes[, c('dataset', 'type')]), all = T)

ggplot(slopes[!(slopes$dataset == 'Bacterial starch hydrolysis' & slopes$species == 'P'), ],
       aes(x = slope_eps, y = empirical_slope, color = dataset, shape = type)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_blank(aes(x = empirical_slope, y = slope_eps)) + # to make x and y scales the same
  geom_point(cex = 3) +
  geom_text(data = mycor,
            aes(x = max_eps, y = min_eps, label = cor_eps),
            color = 'black',
            hjust = 1,
            vjust = 0) +
  scale_color_manual(values = mycolors,
                     breaks = names(mycolors),
                     name = 'Data set') +
  scale_shape_manual(values = c(16, 3),
                     name = 'Type of data set') +
  scale_x_continuous(name = 'FEE slope\nexpected from averaged interactions',
                     breaks = pretty_breaks(n = 2)) +
  scale_y_continuous(name = 'Empirical FEE slope',
                     breaks = pretty_breaks(n = 2)) +
  facet_wrap(~ dataset, nrow = 1,
             scales = 'free') +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = 'none')

ggsave(filename = '../plots/slopes_predictions_eps_bydataset.pdf',
       device = 'pdf',
       dpi = 600,
       width = 300,
       height = 100,
       units = 'mm',
       limitsize = F)

# same but for slopes estimated from phi_ij
lims <- range(c(slopes$empirical_slope, slopes$slope_phi), na.rm = T)
ggplot(slopes[!(slopes$dataset == 'Bacterial starch hydrolysis' & slopes$species == 'P'), ],
       aes(x = slope_phi, y = empirical_slope, color = dataset, shape = type)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_point(cex = 3) +
  scale_color_manual(values = mycolors,
                     breaks = names(mycolors),
                     name = 'Data set') +
  scale_shape_manual(values = c(16, 3),
                     name = 'Type of data set') +
  scale_x_continuous(name = 'FEE slope\nexpected from one- and two-species communities',
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

ggsave(filename = '../plots/slopes_predictions_phi.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

ggplot(slopes[!(slopes$dataset == 'Bacterial starch hydrolysis' & slopes$species == 'P'), ],
       aes(x = slope_phi, y = empirical_slope, color = dataset, shape = type)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = '#d1d3d4') +
  geom_blank(aes(x = empirical_slope, y = slope_phi)) + # to make x and y scales the same
  geom_point(cex = 3) +
  geom_text(data = mycor,
            aes(x = max_phi, y = min_phi, label = cor_phi),
            color = 'black',
            hjust = 1,
            vjust = 0) +
  scale_color_manual(values = mycolors,
                     breaks = names(mycolors),
                     name = 'Data set') +
  scale_shape_manual(values = c(16, 3),
                     name = 'Type of data set') +
  scale_x_continuous(name = 'FEE slope\nexpected from one- and two-species communities',
                     breaks = pretty_breaks(n = 2)) +
  scale_y_continuous(name = 'Empirical FEE slope',
                     breaks = pretty_breaks(n = 2)) +
  facet_wrap(~ dataset, nrow = 1,
             scales = 'free') +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.position = 'none')

ggsave(filename = '../plots/slopes_predictions_phi_bydataset.pdf',
       device = 'pdf',
       dpi = 600,
       width = 300,
       height = 100,
       units = 'mm',
       limitsize = F)

