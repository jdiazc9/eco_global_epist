# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(tidyverse)

# load data
data <- read.csv('../pyoverdine_data/training_set.csv')
data <- lapply(1:3,
               FUN = function(i) matrix2string(data[1:(nrow(data) - 1), c(1:8, 8+i)]))

# full species names
sp_names <- setNames(c('Enterobacter sp.',
                       'Pseudomonas sp. 02',
                       'Klebsiella sp.',
                       'Pseudomonas sp. 03',
                       'Pseudomonas sp. 04',
                       'Raoultella sp.',
                       'Pseudomonas sp. 01',
                       'Pseudomonas sp. 05'),
                     paste('sp', 1:8, sep = '_'))

ge_data <- lapply(data, FUN = makeGEdata)
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
                           levels = c("Pseudomonas sp. 01",
                                      "Pseudomonas sp. 02",
                                      "Pseudomonas sp. 03",
                                      "Pseudomonas sp. 04",
                                      "Pseudomonas sp. 05",
                                      "Enterobacter sp.",
                                      "Raoultella sp.",
                                      "Klebsiella sp.")) # ordered from higher to lower function in monoculture

# assemblage of all contributors & best community
data <- merge(merge(data[[1]], data[[2]], by = 'community'),
              data[[3]], by = 'community')
colnames(data)[2:4] <- paste('fun.rep', 1:3, sep = '')
data$fun.mean <- as.numeric(sapply(1:nrow(data),
                                   FUN = function(i) mean(as.numeric(data[i, 2:4]))))
data$fun.sd <- as.numeric(sapply(1:nrow(data),
                                 FUN = function(i) sd(as.numeric(data[i, 2:4]))))

all_contrib <- data[data$community == 'sp_2,sp_4,sp_5,sp_7,sp_8', ]
best_comm <- data[data$fun.mean > all_contrib$fun.mean, ]

# plot histogram of functions
myplot <- data[, c('community', 'fun.mean')] %>% add_row(community = 'range', fun.mean = range(data$fun.mean) + c(-17, 17)) %>%
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
                     name = 'Community\nfunction [uM]') +
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
       filename = '../plots/pyoverdine_histogram.pdf',
       device = 'pdf',
       dpi = 600,
       width = 75,
       height = 75,
       units = 'mm',
       limitsize = F)

# plot inset of panel E
myplot <-
  ggplot(ge_data[ge_data$knock_in == 'Enterobacter sp.', ],
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
                     name = 'Function of ecological background [uM]') +
  scale_y_continuous(breaks = pretty_breaks(n = 3),
                     name = 'dF [uM]',
                     limits = c(-35, 50)) +
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
       filename = '../plots/stitching_method_inset.pdf',
       device = 'pdf',
       dpi = 600,
       width = 75,
       height = 75,
       units = 'mm',
       limitsize = F)

# predicted function of out of sample communities
test_set <- read.csv('../pyoverdine_data/test_set.csv')
test_set <- merge(matrix2string(test_set[, c(1:9)]), matrix2string(test_set[, c(1:8, 10)]),
                  by = 'community', all = T, suffixes = c('.rep1', '.rep2'))

test_set$fun.mean <- rowMeans(test_set[, 2:3])
test_set$fun.sd <- sapply(1:nrow(test_set),
                          FUN = function(i) as.numeric(sd(test_set[i, 2:3])))

data <- data[, c('community', 'fun.mean')]
colnames(data)[2] <- 'fun'

ge_data <- ge_data[, c('background', 'knock_in', 'background_f.mean', 'd_f.mean'), ]
colnames(ge_data) <- gsub('.mean', '', colnames(ge_data))
ge_data$knock_in <- setNames(names(sp_names), sp_names)[as.character(ge_data$knock_in)]

fits <- makeFEEs(ge_data)
eps <- inferAllResiduals(ge_data)

predicted_f <- predictF_fullClosure(test_set$community, data, eps)

pred_obs <- merge(predicted_f, test_set, by = 'community')
pred_obs <- pred_obs[, c('community', 'fun', 'fun.mean', 'fun.sd')]
colnames(pred_obs)[2] <- 'fun_pred'

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
                  alpha = 0.5) +
    geom_point(shape = 16,
               cex = 2) +
    scale_x_continuous(breaks = pretty_breaks(n = 3),
                       name = 'Predicted F [uM]',
                       limits = range) +
    scale_y_continuous(breaks = pretty_breaks(n = 3),
                       name = ' \nObserved F [uM]',
                       limits = range) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 1,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none') +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

print(myplot)
ggsave(myplot,
       filename = '../plots/pyoverdine_predicted_vs_observed_function.pdf',
       device = 'pdf',
       dpi = 600,
       width = 80,
       height = 80,
       units = 'mm',
       limitsize = F)
