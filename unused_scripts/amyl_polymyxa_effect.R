rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(cowplot)

# load data sets
files <- list.files('../data_sets', full.names = T)
files <- files[grepl('amyl', files)]

data <- read.csv(files)
data <- matrix2string(data)
ge_data <- makeGEdata(data)
ge_data$P_bg <- c('No', 'Yes')[1 + grepl('P', ge_data$background)]
ge_data$poly_group <- ge_data$knock_in
ge_data$poly_group[grepl('T', ge_data$background) & ge_data$knock_in == 'P'] <- 'P.T'

sp_names <- setNames(c('B. cereus', 'B. megaterium', 'B. mojavensis', 'P. polymyxa', 'B. subtilis', 'B. thuringiensis'),
                     c('C', 'E', 'M', 'P', 'S', 'T'))
ge_data$knock_in <- sp_names[ge_data$knock_in]

ggplot(ge_data, aes(x = background_f, y = d_f, group = poly_group)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = 'gray') +
  geom_point(color = 'black',
             shape = 16,
             cex = 3) +
  geom_smooth(method = 'lm',
              formula = y~x,
              se = FALSE,
              color = 'black',
              fullrange = TRUE) +
  scale_x_continuous(name = 'F (background)',
                     breaks = pretty_breaks(n = 3),
                     limits = c(-5, 40)) +
  scale_y_continuous(name = 'dF',
                     breaks = pretty_breaks(n = 3),
                     limits = c(-20, 50)) +
  facet_wrap(~ knock_in,
             nrow = 1) +
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

ge_data$poly_group_full <- paste(ge_data$poly_group, ge_data$P_bg, sep = '-')
ge_data$P_bg <- factor(ge_data$P_bg, levels = c('Yes', 'No'))
ggplot(ge_data, aes(x = background_f, y = d_f, group = poly_group_full, color = P_bg)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = 'gray') +
  geom_point(shape = 16,
             cex = 3) +
  geom_smooth(method = 'lm',
              formula = y~x,
              se = FALSE,
              fullrange = FALSE) +
  scale_x_continuous(name = 'F (background)',
                     breaks = pretty_breaks(n = 3),
                     limits = c(-5, 40)) +
  scale_y_continuous(name = 'dF',
                     breaks = pretty_breaks(n = 3),
                     limits = c(-20, 50)) +
  scale_color_manual(name = 'P. polymyxa in background?',
                     values = c('red', 'black')) +
  facet_wrap(~ knock_in,
             nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 0.6,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
