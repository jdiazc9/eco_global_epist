rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)

# load data sets
data <- lapply(list.files('../data_sets/', full.names = T), FUN = function(file) read.csv(file))

# full species names
sp_names <- vector(mode = 'list', length = length(data))

# amylolytic activity
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

# make plots of dF vs. F_background
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
  ge_data <- makeGEdata(matrix2string(data[[i]]), exclude.single.mut = F)
  
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
       filename = '../plots/fig1.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 500,
       units = 'mm',
       limitsize = F)

