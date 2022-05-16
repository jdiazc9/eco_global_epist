rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(tidyverse)

# load data sets
data <- lapply(list.files('../pyoverdine_data/', full.names = T), FUN = function(file) read.table(file, sep = '\t', header = T))

# communities that were not assembled in replicate 1 will be used to test prediction method ('out of sample' = oos)
which_oos <- data[[2]]$community[!(data[[2]]$community %in% data[[1]]$community)]
oos <- merge(data[[2]][data[[2]]$community %in% which_oos, ],
             data[[3]][data[[3]]$community %in% which_oos, ],
             by = 'community',
             suffixes = c('.rep1', '.rep2'))

data <- lapply(data, FUN = function(x) x[!(x$community %in% which_oos), ])

# full species names
sp_names <- setNames(c('E. mori',
                       'P. putida',
                       'K. grimontii',
                       'P. fulva',
                       'P. parafulva',
                       'R. ornithinolytica',
                       'P. cremoricolorata',
                       'P. savastanoi'),
                     1:8)

ge_data <- lapply(data, FUN = makeGEdata, exclude.single.mut = T)
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

ge_data$knock_in <- sp_names[ge_data$knock_in]
ge_data$knock_in <- factor(ge_data$knock_in,
                           levels = c("K. grimontii",
                                      "E. mori",
                                      "R. ornithinolytica",
                                      "P. savastanoi",
                                      "P. parafulva",
                                      "P. fulva",
                                      "P. putida",
                                      "P. cremoricolorata")) # ordered from lower to higher function in monoculture

# plot panel D
myplot <-
  ggplot(ge_data, aes(x = background_f.mean, xmin = background_f.mean - background_f.sd, xmax = background_f.mean + background_f.sd,
                      y = d_f.mean, ymin = d_f.mean - d_f.sd, ymax = d_f.mean + d_f.sd,
                      group = knock_in)) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_point(color = 'black',
               cex = 1.5,
               shape = 16) +
    geom_errorbar(alpha = 0.25) +
    geom_errorbarh(alpha = 0.25) +
    geom_smooth(method = 'lm',
                formula = y~x,
                color = 'firebrick1',
                se = F,
                fullrange = T) +
    scale_x_continuous(breaks = pretty_breaks(n = 3),
                       name = 'Function of ecological background [a.u.]') +
    scale_y_continuous(breaks = pretty_breaks(n = 3),
                       name = 'dF [a.u.]') +
    facet_wrap(~knock_in,
               nrow = 2) +
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
       filename = '../plots/fig2A.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# assemblage of all contributors
data <- merge(merge(data[[1]], data[[2]], by = 'community'),
              data[[3]], by = 'community')
colnames(data)[2:4] <- paste('fun.rep', 1:3, sep = '')
data$fun.mean <- as.numeric(sapply(1:nrow(data),
                                   FUN = function(i) mean(as.numeric(data[i, 2:4]))))
data$fun.sd <- as.numeric(sapply(1:nrow(data),
                                 FUN = function(i) sd(as.numeric(data[i, 2:4]))))

all_contrib <- data[data$community == '2,4,5,7,8', ]
best_comm <- data[data$fun.mean > all_contrib$fun.mean, ]

# plot histogram of functions (panel C)
myplot <- data[, c('community', 'fun.mean')] %>% add_row(community = 'range', fun.mean = range(data$fun.mean) + c(-0.1, 0.1)) %>%
  ggplot(aes(x = fun.mean)) +
    geom_histogram(aes(y = ..count../sum(..count..)),
                   bins = 30,
                   fill = '#39b54a',
                   alpha = 0.25)
binwidth <- layer_data(myplot) %>% mutate(w=xmax-xmin) %>% pull(w) %>% median
myplot <- myplot +
  stat_bin(aes(y = ..count../sum(..count..)),
                 bins = 30,
                 color = 'black',
                 geom = 'step',
                 position = position_nudge(x = -0.5*binwidth)) +
  coord_cartesian(xlim = range(data$fun.mean[1:(nrow(data)-2)]) + c(-0.05, 0.05)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = 'black') +
  geom_vline(xintercept = all_contrib$fun.mean,
             linetype = 'dashed') +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = 'Community\nfunction [a.u.]') +
  scale_y_continuous(breaks = pretty_breaks(n = 3),
                     name = 'Fraction') +
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
       filename = '../plots/fig2B.pdf',
       device = 'pdf',
       dpi = 600,
       width = 75,
       height = 75,
       units = 'mm',
       limitsize = F)

# plot inset of panel E
myplot <-
  ggplot(ge_data[ge_data$knock_in == 'E. mori', ],
         aes(x = background_f.mean, xmin = background_f.mean - background_f.sd, xmax = background_f.mean + background_f.sd,
             y = d_f.mean, ymin = d_f.mean - d_f.sd, ymax = d_f.mean + d_f.sd,
             group = knock_in)) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 2.5,
             shape = 16) +
  geom_smooth(method = 'lm',
              formula = y~x,
              color = 'firebrick1',
              se = F,
              fullrange = T) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = 'Function of ecological background [a.u.]') +
  scale_y_continuous(breaks = pretty_breaks(n = 3),
                     name = 'dF [a.u.]',
                     limits = c(-0.22, 0.25)) +
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
       filename = '../plots/fig2E.pdf',
       device = 'pdf',
       dpi = 600,
       width = 75,
       height = 75,
       units = 'mm',
       limitsize = F)

# predicted function of out of sample communities
oos$fun.mean <- rowMeans(oos[, 2:3])
oos$fun.sd <- sapply(1:nrow(oos),
                     FUN = function(i) as.numeric(sd(oos[i, 2:3])))

data <- data[, c('community', 'fun.mean')]
colnames(data)[2] <- 'fun'

ge_data <- ge_data[, c('background', 'knock_in', 'background_f.mean', 'd_f.mean'), ]
colnames(ge_data) <- gsub('.mean', '', colnames(ge_data))
ge_data$knock_in <- setNames(names(sp_names), sp_names)[as.character(ge_data$knock_in)]

fits <- makeGEfits(ge_data)
eps <- inferEps(ge_data)

predicted_f <- predictF(oos$community, data, fits, eps)
predicted_f <- data.frame(community = names(predicted_f),
                          fun_pred = as.numeric(predicted_f))

pred_obs <- merge(predicted_f, oos, by = 'community')

r_squared <- cor(pred_obs$fun_pred, pred_obs$fun.mean)^2

# plot predicted vs observed for out of sample communities (panel F)
range <- c(min(c(pred_obs$fun_pred, pred_obs$fun.mean - pred_obs$fun.sd)),
           max(c(pred_obs$fun_pred, pred_obs$fun.mean + pred_obs$fun.sd)))

myplot <- 
  ggplot(pred_obs, aes(x = fun_pred, y = fun.mean)) +
    geom_abline(slope = 1,
                intercept = 0,
                color = '#d1d3d4') +
    geom_errorbar(aes(ymin = fun.mean - fun.sd, ymax = fun.mean + fun.sd),
                  alpha = 0.25) +
    geom_point(shape = 1,
               cex = 2) +
    scale_x_continuous(breaks = pretty_breaks(n = 3),
                       name = 'Predicted F [a.u.]',
                       limits = range) +
    scale_y_continuous(breaks = pretty_breaks(n = 3),
                       name = ' \nObserved F [a.u.]',
                       limits = range) +
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
       filename = '../plots/fig2F.pdf',
       device = 'pdf',
       dpi = 600,
       width = 100,
       height = 100,
       units = 'mm',
       limitsize = F)
