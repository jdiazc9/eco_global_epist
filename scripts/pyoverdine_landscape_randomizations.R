# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)
library(tidyverse)
library(cowplot)

# load data sets
files <- list.files('../pyoverdine_data', full.names = T, recursive = F, pattern = '.csv')
data <- lapply(files, FUN = function(file) read.csv(file))

data <- lapply(data,
               FUN = function(data) cbind(data[, 1:8], fun = rowMeans(data[, 9:ncol(data)])))
data <- do.call(rbind, data)
data <- matrix2string(data)

# species names
sp_names <- setNames(c('Enterobacter sp.',
                       'Pseudomonas sp. 02',
                       'Klebsiella sp.',
                       'Pseudomonas sp. 03',
                       'Pseudomonas sp. 04',
                       'Raoultella sp.',
                       'Pseudomonas sp. 01',
                       'Pseudomonas sp. 05'),
                     paste('sp', 1:8, sep = '_'))

# plot slopes vs intercepts in pyoverdine data
ge_data <- makeGEdata(data)

coords <- lapply(unique(ge_data$knock_in),
                 FUN = function(sp) {
                   
                   lfit <- lm(formula = d_f ~ background_f,
                              data = ge_data[ge_data$knock_in == sp, ])
                   slope <- summary(lfit)$coefficients[2, 1:2]
                   intercept <- summary(lfit)$coefficients[1, 1:2]
                   
                   return(data.frame(species = sp,
                                     slope = as.numeric(slope[1]),
                                     slope.std = as.numeric(slope[2]),
                                     intercept = as.numeric(intercept[1]),
                                     intercept.std = as.numeric(intercept[2])))
                   
                 })
coords <- do.call(rbind, coords)
coords$species <- sp_names[coords$species]
coords$species <- factor(coords$species, levels = c("Pseudomonas sp. 01",
                                                    "Pseudomonas sp. 02",
                                                    "Pseudomonas sp. 03",
                                                    "Pseudomonas sp. 04",
                                                    "Pseudomonas sp. 05",
                                                    "Enterobacter sp.",
                                                    "Raoultella sp.",
                                                    "Klebsiella sp.")) # ordered from higher to lower function in monoculture

ggplot(coords, aes(x = intercept, y = slope,
                   xmin = intercept - intercept.std, xmax = intercept + intercept.std,
                   ymin = slope - slope.std, ymax = slope + slope.std,
                   color = species)) +
  geom_errorbar(width = 0) +
  geom_errorbarh(height = 0) +
  geom_point(shape = 16,
             size = 3) +
  scale_color_manual(values = c('#8b8131',
                                '#ddc85d',
                                '#518b3f',
                                '#86c65e',
                                '#5ea5be',
                                '#b199c1',
                                '#c25ea4',
                                '#662b85',
                                'black'),
                     name = '') +
  scale_x_continuous(name = 'FEE intercept [uM]',
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = 'FEE slope',
                     breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/pyoverdine_fits_coefficients.pdf',
       device = 'pdf',
       dpi = 600,
       width = 160,
       height = 80,
       units = 'mm',
       limitsize = F)


### LANDSCAPE RANDOMIZATIONS

# wrapper function: randomize data, return FEE coefficients for each species and r_squared of predicted vs. observed function
randomFEEs <- function(data) {
  
  # FEEs in randomized data
  data.rnd <- data.frame(community = data$community,
                         fun = sample(data$fun, size = nrow(data), replace = F))
  ge_data <- makeGEdata(data.rnd)
  fits <- makeFEEs(ge_data)
  
  fits <- data.frame(species = rownames(fits),
                     a = fits[, 'a'],
                     b = fits[, 'b'])
  
  # leave 20% of data out of sample
  which_oos <- sample(1:nrow(data.rnd), size = round(0.2*nrow(data.rnd)), replace = F)
  oos.rnd <- data.rnd[which_oos, ]
  data.rnd <- data.rnd[-which_oos, ]
  oos <- data[which_oos, ]
  data <- data[-which_oos, ]
  if ('' %in% oos$community) {
    data.rnd <- rbind(data.rnd, oos.rnd[oos.rnd$community == '', ])
    data <- rbind(data, oos[oos$community == '', ])
    oos.rnd <- oos[oos.rnd$community != '', ]
    oos <- oos[oos$community != '', ]
  }
  
  # predict using randomized data
  ge_data.rnd <- makeGEdata(data.rnd)
  eps.rnd <- inferAllResiduals(ge_data.rnd)
  predF.rnd <- predictF_fullClosure(oos.rnd$community, data.rnd, eps.rnd)
  
  po.rnd <- merge(predF.rnd, oos.rnd, by = 'community', suffixes = c('_pred', '_obs'))
  r_squared.rnd <- cor(po.rnd$fun_pred, po.rnd$fun_obs)^2
  
  # predict using true data
  ge_data <- makeGEdata(data)
  eps <- inferAllResiduals(ge_data)
  predF <- predictF_fullClosure(oos$community, data, eps)
  
  po <- merge(predF, oos, by = 'community', suffixes = c('_pred', '_obs'))
  r_squared <- cor(po$fun_pred, po$fun_obs)^2
  
  # output
  return(list(fits = fits,
              po = po,
              po.rnd = po.rnd,
              r_squared = r_squared,
              r_squared.rnd = r_squared.rnd))
  
}

# do 500 randomizations
rsq <- NULL
rsq.rnd <- NULL
fits <- data.frame(run = numeric(0),
                   species = character(0),
                   a = numeric(0),
                   b = numeric(0))
for (i in 1:500) {
  
  print(i)
  
  run <- randomFEEs(data)
  rsq <- c(rsq, run$r_squared)
  rsq.rnd <- c(rsq.rnd, run$r_squared.rnd)
  fits <- rbind(fits, cbind(run = i, run$fits))
  
}

# plot FEEs in randomized data
data.rnd <- data.frame(community = data$community,
                       fun = sample(data$fun, size = nrow(data), replace = F))
ge_data.rnd <- makeGEdata(data.rnd)
ge_data.rnd$knock_in <- sp_names[ge_data.rnd$knock_in]
ge_data.rnd$knock_in <- factor(ge_data.rnd$knock_in, levels = c("Pseudomonas sp. 01",
                                                                "Pseudomonas sp. 02",
                                                                "Pseudomonas sp. 03",
                                                                "Pseudomonas sp. 04",
                                                                "Pseudomonas sp. 05",
                                                                "Enterobacter sp.",
                                                                "Raoultella sp.",
                                                                "Klebsiella sp.")) # ordered from higher to lower function in monoculture


ggplot(ge_data.rnd, aes(x = background_f, y = d_f)) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             cex = 1.5,
             shape = 16) +
  geom_smooth(method = 'lm',
              formula = y~x,
              color = 'firebrick1',
              se = F,
              fullrange = T) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     name = 'Function of ecological background [uM]') +
  scale_y_continuous(breaks = pretty_breaks(n = 3),
                     name = 'dF [uM]') +
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

ggsave(filename = '../plots/pyoverdine_landscape_randomizations_FEEs.pdf',
       device = 'pdf',
       dpi = 600,
       width = 200,
       height = 100,
       units = 'mm',
       limitsize = F)

# plot a vs. b (true vs randomized landscapes)
fits.true <- makeFEEs(makeGEdata(data))
fits <- rbind(data.frame(run = 'true',
                         species = rownames(fits.true),
                         a = fits.true[, 'a'],
                         b = fits.true[, 'b']),
              fits)
fits$color <- fits$species
fits$color[fits$run != 'true'] <- 'Randomized data'
fits$color[fits$run == 'true'] <- sp_names[fits$color[fits$run == 'true']]
fits$color <- factor(fits$color, levels = c("Pseudomonas sp. 01",
                                            "Pseudomonas sp. 02",
                                            "Pseudomonas sp. 03",
                                            "Pseudomonas sp. 04",
                                            "Pseudomonas sp. 05",
                                            "Enterobacter sp.",
                                            "Raoultella sp.",
                                            "Klebsiella sp.",
                                            "Randomized data"))

ggplot(fits, aes(x = a, y = b, color = color, alpha = color, size = color)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c('#8b8131',
                                '#ddc85d',
                                '#518b3f',
                                '#86c65e',
                                '#5ea5be',
                                '#b199c1',
                                '#c25ea4',
                                '#662b85',
                                'black'),
                     name = '') +
  scale_alpha_manual(values = c(rep(1, 8), 0.2),
                     guide = 'none') +
  scale_size_manual(values = c(rep(3, 8), 2),
                    guide = 'none') +
  scale_x_continuous(name = 'FEE intercept [uM]',
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = 'FEE slope',
                     breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic',
                                  size = 10),
        aspect.ratio = 1,
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/pyoverdine_landscape_randomizations_fits_coefficients.pdf',
       device = 'pdf',
       dpi = 600,
       width = 160,
       height = 80,
       units = 'mm',
       limitsize = F)

# plot histograms of R^2 pred vs. obs
r <- data.frame(run = c(rep('True data', length(rsq)), rep('Randomized data', length(rsq.rnd))),
                rsq = c(rsq, rsq.rnd))
r$run <- factor(r$run, levels = c('True data', 'Randomized data'))

myplot <-
  ggplot(r, aes(x = rsq, fill = run)) +
    geom_histogram(bins = 50,
                   color = 'black') +
    scale_fill_manual(values = c('white', 'black'),
                      name = '') +
    scale_x_continuous(name = expression(paste(italic(R)^2, ' predicted vs. observed')),
                       breaks = c(0, 0.5, 1),
                       labels = c('0', '0.5', '1'),
                       limits = c(-0.1, 1.1)) +
    scale_y_continuous(name = 'Frequency',
                       expand = expansion(mult = c(0, .1))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'italic',
                                    size = 10),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 14)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

print(myplot)
ggsave(myplot,
       filename = '../plots/pyoverdine_landscape_randomizations_predictions_histograms.pdf',
       device = 'pdf',
       dpi = 600,
       width = 160,
       height = 80,
       units = 'mm',
       limitsize = F)
