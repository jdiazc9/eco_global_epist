rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)

# load amylase data set
data <- read.csv('../data_sets/amyl_Sanchez-Gorostiaga2019.csv')

# full species names
sp_names <- setNames(c('B. cereus', 'B. megaterium', 'B. mojavensis', 'P. polymyxa', 'B. subtilis', 'B. thuringiensis'),
                     c('C', 'E', 'M', 'P', 'S', 'T'))

# GE data
data <- makeGEdata(matrix2string(data))

# focus on P. polymyxa
data <- data[data$knock_in == 'P', ]

# is B. thuringiensis in the ecological background?
data$B_background <- sapply(data$background, FUN = function(x) grepl('T', x))
data$B_background <- c('No', 'Yes')[1 + data$B_background]
data$B_background <- factor(data$B_background, levels = c('Yes', 'No'))

# plot
myplot <-
  ggplot(data, aes(x = background_f, y = d_f, shape = B_background, group = B_background)) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_point(cex = 3) +
    geom_smooth(method = 'lm',
                se = F,
                fullrange = T,
                formula = y~x,
                color = 'deepskyblue') +
    scale_x_continuous(name = 'Function of ecological background [a.u.]',
                       breaks = pretty_breaks(n = 3),
                       limits = c(-5, 20)) +
    scale_y_continuous(name = 'dF [a.u.]',
                       breaks = pretty_breaks(n = 3),
                       limits = c(-10, 40)) +
    scale_shape_manual(values = c(16, 1),
                       name = expression(paste(italic('B. thuringiensis'), ' in background?'))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    ggtitle(expression(italic('P. polymyxa'))) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

print(myplot)
ggsave(myplot,
       filename = '../plots/figS1.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 200,
       units = 'mm',
       limitsize = F)

