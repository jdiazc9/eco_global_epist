# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)

# load data sets
files <- list.files('../data_sets', full.names = T)
files <- c(files, '../pyoverdine_data/training_set.csv')
data <- lapply(files, FUN = function(file) read.csv(file))

data[[6]][, 9] <- rowMeans(data[[6]][, c(9, 10, 11)])
data[[6]]$function.rep2 <- NULL
data[[6]]$function.rep3 <- NULL
colnames(data[[6]])[9] <- 'function.'

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


# monocultures
f_monoc <- lapply(1:length(data),
                  FUN = function(i) {
                    
                    datai <- matrix2string(data[[i]])
                    datai <- aggregate(fun ~ community,
                                       data = datai,
                                       FUN = mean)
                    
                    f_monoc <- datai[nSpecies(datai$community) == 1, ]
                    f_monoc$community <- sp_names[[i]][f_monoc$community]
                    return(f_monoc)
                    
                  })

# functional optima
f_opt <- lapply(1:length(data),
                FUN = function(i) {
                  
                  datai <- data[[i]]
                  colnames(datai)[1:(ncol(datai) - 1)] <- sp_names[[i]][colnames(datai)[1:(ncol(datai) - 1)]]
                  datai <- matrix2string(datai)
                  datai <- aggregate(fun ~ community,
                                     data = datai,
                                     FUN = mean)
                  
                  f_opt <- datai[which.max(datai$fun), ]
                  f_opt$community <- gsub(',', '\n', f_opt$community)
                  return(f_opt)
                  
                })

# scale monoculture functions by optimal functions and add data set names
data_set <- basename(files)

for (i in 1:length(f_monoc)) {
  
  f_monoc[[i]]$fun <- f_monoc[[i]]$fun/f_opt[[i]]$fun
  f_opt[[i]]$fun <- 1
  
  f_monoc[[i]] <- cbind(data_set = data_set[i],
                        is_monoc = 'yes',
                        color = 'monoc',
                        f_monoc[[i]])
  f_opt[[i]] <- cbind(data_set = data_set[i],
                      is_monoc = c('yes', 'no')[1 + grepl('\\n', f_opt[[i]]$community)],
                      color = 'opt',
                      f_opt[[i]])
  
}

# bind data sets
f <- rbind(do.call(rbind, f_opt),
           do.call(rbind, f_monoc))
remove_these <- which(f$color == 'opt' & f$is_monoc == 'yes')
f <- f[-remove_these, ]

f$data_set <- setNames(c('Bacterial\nstarch hydrolysis',
                         'Bacterial\nbutyrate secretion',
                         'Phytoplankton\nbiomass',
                         'Above-ground\nplant biomass',
                         'Bacterial\nxylose oxidation',
                         'Bacterial\npyoverdines secretion'),
                       basename(files))[f$data_set]

f$data_set <- factor(f$data_set, levels = c('Bacterial\npyoverdines secretion',
                                            'Above-ground\nplant biomass',
                                            'Phytoplankton\nbiomass',
                                            'Bacterial\nxylose oxidation',
                                            'Bacterial\nstarch hydrolysis',
                                            'Bacterial\nbutyrate secretion'))
f$community[f$color == 'monoc'] <- ''

# plot
ggplot(f, aes(x = 0, y = fun, color = color)) +
  geom_hline(yintercept = 1,
             color = '#d1d3d4') +
  geom_jitter(width = 0.25,
              size = 3) +
  geom_text(aes(label = community),
            fontface = 'italic',
            size = 5,
            hjust = 0,
            nudge_y = -0.25,
            nudge_x = -0.5) +
  facet_wrap(~data_set,
             nrow = 1) +
  scale_x_continuous(name = '',
                     limits = c(-1, 1)) +
  scale_y_continuous(name = expression(paste(italic(F), ' / ', italic(F)[max])),
                     limits = c(-0.01, 1.01),
                     breaks = c(0, 0.5, 1),
                     labels = c('0', '0.5', '1')) +
  scale_color_manual(values = c('black', 'red'),
                     name = '',
                     labels = c('Monocultures',
                                expression(paste(n, ' \u2265 ', '2 species')))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        aspect.ratio = 3,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5) +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf,size=0.5) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/datasets_functional_maxima.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)





