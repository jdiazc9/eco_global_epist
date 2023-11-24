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
ge_data <- lapply(ge_data,
                  FUN = function(df) {
                    df$knockin_f <- df$background_f + df$d_f
                    return(df)
                  })
ge_data <- merge(merge(ge_data[[1]], ge_data[[2]],
                       by = c('background', 'knock_in'),
                       suffixes = c('.rep1', '.rep2'),
                       all = T),
                 ge_data[[3]],
                 by = c('background', 'knock_in'),
                 suffixes = c('', '.rep3'),
                 all = T)
colnames(ge_data)[9:11] <- c('background_f.rep3', 'd_f.rep3', 'knockin_f.rep3')

ge_data$background_f.mean <- as.numeric(sapply(1:nrow(ge_data),
                                               FUN = function(i) mean(as.numeric(ge_data[i, paste('background_f.rep', 1:3, sep = '')]))))
ge_data$background_f.sd <- as.numeric(sapply(1:nrow(ge_data),
                                             FUN = function(i) sd(as.numeric(ge_data[i, paste('background_f.rep', 1:3, sep = '')]))))
ge_data$d_f.mean <- as.numeric(sapply(1:nrow(ge_data),
                                      FUN = function(i) mean(as.numeric(ge_data[i, paste('d_f.rep', 1:3, sep = '')]))))
ge_data$d_f.sd <- as.numeric(sapply(1:nrow(ge_data),
                                    FUN = function(i) sd(as.numeric(ge_data[i, paste('d_f.rep', 1:3, sep = '')]))))
ge_data$knockin_f.mean <- as.numeric(sapply(1:nrow(ge_data),
                                               FUN = function(i) mean(as.numeric(ge_data[i, paste('knockin_f.rep', 1:3, sep = '')]))))
ge_data$knockin_f.sd <- as.numeric(sapply(1:nrow(ge_data),
                                             FUN = function(i) sd(as.numeric(ge_data[i, paste('knockin_f.rep', 1:3, sep = '')]))))

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
ge_data_pyo <- ge_data

# fit via total least squares regression
fees_tls <- do.call(rbind,
                    lapply(unique(ge_data_pyo$knock_in),
                           FUN = function(sp) {
                             
                             mytls <- prcomp(ge_data_pyo[ge_data_pyo$knock_in == sp, c('background_f.mean', 'knockin_f.mean')])$rotation
                             slope <- mytls[2, 1]/mytls[1, 1]
                             intercept <- mean(ge_data_pyo$knockin_f.mean[ge_data_pyo$knock_in == sp]) - slope*mean(ge_data_pyo$background_f.mean[ge_data_pyo$knock_in == sp])
                             
                             return(data.frame(knock_in = sp,
                                               slope = slope,
                                               intercept = intercept))
                             
                           }))

# plot FEEs (in F vs. F format)
ggplot(ge_data_pyo, aes(x = background_f.mean, xmin = background_f.mean - background_f.sd, xmax = background_f.mean + background_f.sd,
                    y = knockin_f.mean, ymin = knockin_f.mean - knockin_f.sd, ymax = knockin_f.mean + knockin_f.sd,
                    group = knock_in)) +
  geom_abline(slope = 1, intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 2,
             shape = 1) +
  geom_errorbar(alpha = 0.25) +
  geom_errorbarh(alpha = 0.25) +
  # geom_smooth(method = 'lm',
  #             formula = y~x,
  #             color = 'firebrick1',
  #             se = F,
  #             fullrange = F) +
  geom_abline(data = fees_tls,
              aes(slope = slope,
                  intercept = intercept,
                  color = slope), size = 1) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = 'Function of ecological background [uM]') +
  scale_y_continuous(breaks = pretty_breaks(n = 3),
                     name = 'Function with focal species [uM]') +
  scale_color_gradientn(colors = c('firebrick1',
                                   'firebrick1',
                                   'black',
                                   'deepskyblue',
                                   'deepskyblue'),
                        breaks = 1 + c(-0.8, -0.5, 0, 0.5, 0.8),
                        limits = 1 + c(-0.8, 0.8),
                        na.value = 'deepskyblue') +
  facet_wrap(~knock_in,
             nrow = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 16),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/fees_pyoverdine_FvsF.pdf',
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
ge_data$knockin_f <- ge_data$background_f + ge_data$d_f
ge_data_khan <- ge_data

fees_tls <- do.call(rbind,
                    lapply(unique(ge_data_khan$knock_in),
                           FUN = function(sp) {
                             
                             mytls <- prcomp(ge_data_khan[ge_data_khan$knock_in == sp, c('background_f', 'knockin_f')])$rotation
                             slope <- mytls[2, 1]/mytls[1, 1]
                             intercept <- mean(ge_data_khan$knockin_f[ge_data_khan$knock_in == sp]) - slope*mean(ge_data_khan$background_f[ge_data_khan$knock_in == sp])
                             
                             return(data.frame(knock_in = sp,
                                               slope = slope,
                                               intercept = intercept))
                             
                           }))

ggplot(ge_data_khan, aes(x = background_f, y = knockin_f, group = knock_in)) +
  geom_abline(slope = 1, intercept = 0,
              color = '#d1d3d4') +
  geom_blank(aes(x = knockin_f, y = background_f)) + # to homogenize x and y scales
  geom_point(color = 'black',
             cex = 3,
             shape = 1) +
  # geom_smooth(method = 'lm',
  #             formula = y~x,
  #             color = 'firebrick1',
  #             se = F,
  #             fullrange = F) +
  geom_abline(data = fees_tls,
              aes(slope = slope,
                  intercept = intercept,
                  color = slope), size = 1) +
  scale_color_gradientn(colors = c('firebrick1',
                                   'firebrick1',
                                   'black',
                                   'deepskyblue',
                                   'deepskyblue'),
                        breaks = 1 + c(-0.8, -0.5, 0, 0.5, 0.8),
                        limits = 1 + c(-0.8, 0.8),
                        na.value = 'deepskyblue') +
  scale_x_continuous(breaks = pretty_breaks(n = 2),
                     name = expression(paste('Relative fitness of genetic background, ', italic(f)[B], sep = ''))) +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = 'Fitness with\nfocal mutation') +
  facet_wrap(~knock_in,
             nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 16),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/fees_khan_FvsF.pdf',
       device = 'pdf',
       dpi = 600,
       width = 160,
       height = 100,
       units = 'mm',
       limitsize = F)

# make plots of F_knockin vs. F_background (all other datasets, this goes in fig. 2)
g <- vector(mode = 'list', length = length(data))

for (i in 1:5) {
  
  # if there are multiple measurements of a same community, average them
  colnames(data[[i]])[ncol(data[[i]])] <- 'fun'
  data[[i]] <- aggregate(fun ~ .,
                         data = data[[i]],
                         FUN = mean)
  
  # for the phytoplankton biomass dataset (Ghedini et al., scale functions by 1e-4 for easier readability)
  if (i == 3) data[[i]][, ncol(data[[i]])] <- data[[i]][, ncol(data[[i]])]/1e4
  
  # extract F vs. dF data from combinatorial assemblages
  ge_data <- makeGEdata(matrix2string(data[[i]]))
  ge_data$knockin_f <- ge_data$background_f + ge_data$d_f
  
  # remove species for which there are not at least 5 data points (most data sets are combinatorially incomplete)
  npoints <- sapply(unique(ge_data$knock_in),
                    FUN = function(sp) sum(ge_data$knock_in == sp))
  valid_sp <- names(npoints)[npoints >= 5]
  ge_data <- ge_data[ge_data$knock_in %in% valid_sp, ]
  
  # group for linear fits (group by species unless branching is observed, i.e. P. polymyxa when B. thruingiensis is present in amyl data)
  ge_data$group <- ge_data$knock_in
  
  if (any(ge_data$knock_in == 'P')) ge_data$group[ge_data$knock_in == 'P' & grepl('T', ge_data$background)] <- 'P.T'
  
  # make total-least squares fits and get slopes
  fees_tls <- do.call(rbind,
                      lapply(unique(ge_data$group),
                             FUN = function(sp) {
                               
                               mytls <- prcomp(ge_data[ge_data$group == sp, c('background_f', 'knockin_f')])$rotation
                               slope <- mytls[2, 1]/mytls[1, 1]
                               intercept <- mean(ge_data$knockin_f[ge_data$group == sp]) - slope*mean(ge_data$background_f[ge_data$group == sp])
                               
                               return(data.frame(knock_in = sp,
                                                 slope = slope,
                                                 intercept = intercept))
                               
                             }))
  fees_tls$color <- fees_tls$slope
  fees_tls$color[fees_tls$color > 1.8] <- 1.8 # 'shrink' color values from the tails of the distribution to avoid NAs in color scale
  fees_tls$color[fees_tls$color < 0.2] <- 0.2
  
  # manually set axis limits for clear visualization
  dx <- (max(ge_data$background_f) - min(ge_data$background_f))/20
  dy <- (max(ge_data$knockin_f) - min(ge_data$knockin_f))/20
  
  dx <- c(min(ge_data$background_f) - dx, max(ge_data$background_f) + dx)
  dy <- c(min(ge_data$knockin_f) - dy, max(ge_data$knockin_f) + dy)
  
  # full species names
  ge_data$knock_in <- sp_names[[i]][ge_data$knock_in]
  fees_tls$knock_in <- setNames(unique(ge_data[, c('knock_in', 'group')])$knock_in,
                                unique(ge_data[, c('knock_in', 'group')])$group)[fees_tls$knock_in]
  
  # if there are less than 25 species (maximum across our data sets), add 'ghost' panels to mantain axis proportions
  if (length(unique(ge_data$knock_in)) < 25) {
    
    add_panels <- 25 - length(unique(ge_data$knock_in))
    
    ge_data <- rbind(ge_data,
                     data.frame(background = NA,
                                knock_in = paste('remove_this_panel', 1:add_panels, sep = '.'),
                                background_f = NA,
                                knockin_f = NA,
                                d_f = NA,
                                group = NA))
    
  }
  ge_data$knock_in <- factor(ge_data$knock_in, levels = unique(ge_data$knock_in))
  fees_tls$knock_in <- factor(fees_tls$knock_in, levels = levels(ge_data$knock_in))
  
  # plot
  g[[i]] <- 
    ggplot(ge_data, aes(x = background_f, y = knockin_f, group = group)) +
    geom_abline(slope = 1,
                intercept = 0,
                color = 'gray') +
    geom_blank(aes(x = knockin_f, y = background_f)) + # to homogenize y and x scales
    geom_point(color = 'black',
               shape = 1,
               cex = 2) +
    # geom_smooth(method = 'lm',
    #             formula = y~x,
    #             se = FALSE,
    #             fullrange = FALSE) +
    geom_abline(data = fees_tls,
                aes(slope = slope,
                    intercept = intercept,
                    color = color), size = 1) +
    scale_x_continuous(name = '',
                       breaks = pretty_breaks(n = 2)) +
    scale_y_continuous(name = 'F (knockin)',
                       breaks = pretty_breaks(n = 2)) +
    facet_wrap(~ knock_in,
               ncol = 7) +
    scale_color_gradientn(colors = c('firebrick1',
                                     'firebrick1',
                                     'black',
                                     'deepskyblue',
                                     'deepskyblue'),
                          breaks = 1 + c(-0.8, -0.5, 0, 0.5, 0.8),
                          limits = 1 + c(-0.8, 0.8),
                          na.value = 'yellow') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 1,
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
       filename = '../plots/datasets_fees_FvsF.pdf',
       device = 'pdf',
       dpi = 600,
       width = 300,
       height = 600,
       units = 'mm',
       limitsize = F)

# quality of fits
# for species i, quality is quantified as the R2 of the F(knockin) vs. F(background) plot

sigma <- data.frame(dataset = character(0),
                    species = character(0),
                    R2_FvsF = numeric(0))

for (i in 1:6) {
  
  ge_data <- makeGEdata(matrix2string(data[[i]]))
  
  ge_data$group <- ge_data$knock_in
  if (any(ge_data$knock_in == 'P')) ge_data$group[ge_data$knock_in == 'P' & grepl('T', ge_data$background)] <- 'P.T'
  ge_data$knockin_f <- ge_data$background_f + ge_data$d_f
  
  # make linear fits
  lfits <- sapply(unique(ge_data$group),
                  FUN = function(grp) summary(lm(formula = knockin_f ~ background_f,
                                                 data = ge_data[ge_data$group == grp, ]))$r.squared)
  fees_tls <- do.call(rbind,
                      lapply(unique(ge_data$group),
                             FUN = function(sp) {
                               
                               mytls <- prcomp(ge_data[ge_data$group == sp, c('background_f', 'knockin_f')])$rotation
                               slope <- mytls[2, 1]/mytls[1, 1]
                               intercept <- mean(ge_data$knockin_f[ge_data$group == sp]) - slope*mean(ge_data$background_f[ge_data$group == sp])
                               
                               x <- ge_data$background_f[ge_data$group == sp]
                               y <- ge_data$knockin_f[ge_data$group == sp]
                               yfit <- intercept + slope*x
                               R2 <- 1 - sum((y - yfit)^2)/sum((y - mean(y))^2)
                               R2[R2 < 0] <- 0 # if total-least-squares regression model fails, its R2 is set to zero
                               
                               return(data.frame(knock_in = sp,
                                                 slope = slope,
                                                 intercept = intercept,
                                                 R2 = R2))
                               
                             }))
  
  sigma <- rbind(sigma,
                 data.frame(dataset = basename(files)[i],
                            species = fees_tls$knock_in,
                            R2_FvsF = fees_tls$R2))
  
}

sigma$dataset <- factor(sigma$dataset,
                        levels = c("training_set.csv",
                                   "plant-biommass_Kuebbing2016_all.csv",
                                   "phytoplankton-biomass_Ghedini2022.csv",
                                   "xylose_Langenheder2010.csv",
                                   "amyl_Sanchez-Gorostiaga2019.csv",
                                   "butyrate_Clark2021.csv"))

# plot
myplot1 <- 
  ggplot(sigma, aes(x = 0, y = R2_FvsF)) +
  geom_jitter(size = 3) +
  facet_wrap(~dataset, nrow = 1) +
  scale_x_continuous(name = '',
                     limits = c(-0.75, 0.75)) +
  scale_y_continuous(name = expression(paste(italic(R)^2, ' of ', italic(F)(italic(B)+italic(i)), ' vs. ', italic(F)(italic(B)), sep = '')),
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
  ggplot(sigma, aes(x = R2_FvsF)) +
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
       filename = '../plots/datasets_fees_quality_FvsF.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# save FEE parameters fit to the FvsF representation
fees_fvsf <- do.call(rbind,
                     lapply(1:6,
                            FUN = function(i) {
                              
                              ge_data <- makeGEdata(matrix2string(data[[i]]))
                              
                              ge_data$group <- ge_data$knock_in
                              if (any(ge_data$knock_in == 'P')) ge_data$group[ge_data$knock_in == 'P' & grepl('T', ge_data$background)] <- 'P.T'
                              ge_data$knockin_f <- ge_data$background_f + ge_data$d_f
                              
                              # make total least squares fits
                              fees_tls <- do.call(rbind,
                                                  lapply(unique(ge_data$group),
                                                         FUN = function(sp) {
                                                           
                                                           mytls <- prcomp(ge_data[ge_data$group == sp, c('background_f', 'knockin_f')])$rotation
                                                           slope <- mytls[2, 1]/mytls[1, 1]
                                                           intercept <- mean(ge_data$knockin_f[ge_data$group == sp]) - slope*mean(ge_data$background_f[ge_data$group == sp])
                                                           
                                                           x <- ge_data$background_f[ge_data$group == sp]
                                                           y <- ge_data$knockin_f[ge_data$group == sp]
                                                           yfit <- intercept + slope*x
                                                           R2 <- 1 - sum((y - yfit)^2)/sum((y - mean(y))^2)
                                                           R2[R2 < 0] <- 0 # if total-least-squares regression model fails, its R2 is set to zero
                                                           
                                                           return(data.frame(dataset = basename(files)[i],
                                                                             species = sp,
                                                                             slope = slope,
                                                                             intercept = intercept,
                                                                             R2 = R2))
                                                           
                                                         }))
                              
                              if(grepl('amyl', basename(files)[i])) {
                                fees_tls$species[fees_tls$species == 'P'] <- 'P_0'
                                fees_tls$species[fees_tls$species == 'P.T'] <- 'P_1'
                              }
                              if(grepl('training_set', basename(files)[i])) fees_tls$dataset <- 'pyo'
                              
                              return(fees_tls)
                              
                            }))

# check if slopes differ from one
fees_fvsf$dataset <- factor(fees_fvsf$dataset,
                            levels = c("pyo",
                                       "plant-biommass_Kuebbing2016_all.csv",
                                       "phytoplankton-biomass_Ghedini2022.csv",
                                       "xylose_Langenheder2010.csv",
                                       "amyl_Sanchez-Gorostiaga2019.csv",
                                       "butyrate_Clark2021.csv"))

ggplot(fees_fvsf, aes(x = 0, y = slope)) +
  geom_abline(slope = 0, intercept = 1, color = 'gray') +
  geom_blank(aes(x = 0, y = 1)) +
  geom_point(cex = 2,
             position = position_jitter(width = 0.5)) +
  facet_wrap(~dataset, nrow = 1, scales = 'free_y') +
  scale_x_continuous(name = '',
                     limits = c(-0.75, 0.75)) +
  scale_y_continuous(name = expression(paste('Slope of ', italic(F)(bold(s)+bold(i)), ' vs. ', italic(F)(bold(s)), ' regression', sep = ''))) +
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

ggsave(filename = '../plots/datasets_slopes_FvsF.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# compare with FEEs fit to the dF-vs-F representation
fees <- do.call(rbind,
                lapply(1:6,
                       FUN = function(i) {
                         
                         ge_data <- makeGEdata(matrix2string(data[[i]]))
                         
                         ge_data$group <- ge_data$knock_in
                         if (any(ge_data$knock_in == 'P')) ge_data$group[ge_data$knock_in == 'P' & grepl('T', ge_data$background)] <- 'P.T'
                         ge_data$knockin_f <- ge_data$background_f + ge_data$d_f
                         
                         # make linear fits
                         lfits <- do.call(rbind,
                                          lapply(unique(ge_data$group),
                                                 FUN = function(grp) {
                                                   lfit_i <- summary(lm(formula = d_f ~ background_f,
                                                                        data = ge_data[ge_data$group == grp, ]))
                                                   out <- data.frame(dataset = basename(files)[i],
                                                                     species = grp,
                                                                     slope = lfit_i$coefficients[2, 1],
                                                                     intercept = lfit_i$coefficients[1, 1],
                                                                     slope_stderr = lfit_i$coefficients[2, 2],
                                                                     intercept_stderr = lfit_i$coefficients[1, 2])
                                                   return(out)
                                                 }))
                         
                         if(grepl('amyl', basename(files)[i])) {
                           lfits$species[lfits$species == 'P'] <- 'P_0'
                           lfits$species[lfits$species == 'P.T'] <- 'P_1'
                         }
                         if(grepl('training_set', basename(files)[i])) lfits$dataset <- 'pyo'
                         
                         return(lfits)
                         
                       }))

slopes_cmp <- merge(fees[, c('dataset', 'species', 'slope')],
                    fees_fvsf[, c('dataset', 'species', 'slope')],
                    by = c('dataset', 'species'),
                    suffixes = c('', '_fvsf'))
intercept_cmp <- merge(fees[, c('dataset', 'species', 'intercept')],
                       fees_fvsf[, c('dataset', 'species', 'intercept')],
                       by = c('dataset', 'species'),
                       suffixes = c('', '_fvsf'))

ggplot(slopes_cmp, aes(x = slope, y = slope_fvsf)) +
  geom_abline(slope = 1, 
              intercept = 1,
              color = 'gray') +
  geom_point() +
  geom_blank(aes(x = slope_fvsf, y = slope)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

ggplot(intercept_cmp, aes(x = intercept, y = intercept_fvsf)) +
  geom_point()







