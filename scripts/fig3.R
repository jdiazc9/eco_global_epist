rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)

data <- read.csv('../crm_data/CRM_Succinate_final.txt')
data <- matrix2string(data)
data$fun <- data$fun*100 # re-scale for easier visualization
ge_data <- makeGEdata(data)

ge_data$knock_in <- gsub('S', 'Sp. ', ge_data$knock_in)
ge_data$knock_in <- factor(ge_data$knock_in, levels = paste('Sp.', 0:67))

# reduce size of ge_data do that figure S7 doesn't end up being too heavy
set.seed(0)
ge_data_small <- do.call(rbind,
                         lapply(unique(ge_data$knock_in),
                                FUN = function(sp) {
                                  ge_data_i <- ge_data[ge_data$knock_in == sp, ]
                                  ge_data_i <- ge_data_i[sample(1:nrow(ge_data_i), replace = F, size = 200), ]
                                  return(ge_data_i)
                                }))

# full plot
myplot <-
  ggplot(ge_data_small, aes(x = background_f, y = d_f)) +
    facet_wrap(~knock_in,
               ncol = 7) +
    geom_abline(slope = 0, intercept = 0,
                color = '#d1d3d4') +
    geom_point(color = 'black',
               shape = 1,
               cex = 3) +
    geom_smooth(method = 'lm',
                formula = y~x,
                se = FALSE,
                fullrange = TRUE,
                color = 'firebrick1') +
    scale_x_continuous(name = expression(paste(italic(F), ' (background) [a.u.]')),
                       breaks = pretty_breaks(n = 2)) +
    scale_y_continuous(name = 'dF [a.u.]',
                       breaks = pretty_breaks(n = 2)) +
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

print(myplot)
ggsave(myplot,
       filename = '../plots/figS7.pdf',
       device = 'pdf',
       dpi = 600,
       width = 400,
       height = 400,
       units = 'mm',
       limitsize = F)

# plot of species subset
sp <- paste('Sp.', c(2, 34, 43, 44, 45, 46, 56, 62, 64))

myplot <-
  ggplot(ge_data_small[ge_data_small$knock_in %in% sp, ], aes(x = background_f, y = d_f)) +
  facet_wrap(~knock_in,
             nrow = 3) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_point(color = 'black',
             shape = 1,
             cex = 3,
             alpha = 1) +
  geom_smooth(method = 'lm',
              formula = y~x,
              se = FALSE,
              fullrange = TRUE,
              color = 'firebrick1') +
  scale_x_continuous(name = expression(paste(italic(F), ' (background) [a.u.]')),
                     breaks = pretty_breaks(n = 2)) +
  scale_y_continuous(name = 'dF [a.u.]',
                     breaks = pretty_breaks(n = 2)) +
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

print(myplot)
ggsave(myplot,
       filename = '../plots/fig3A.pdf',
       device = 'pdf',
       dpi = 600,
       width = 160,
       height = 160,
       units = 'mm',
       limitsize = F)


