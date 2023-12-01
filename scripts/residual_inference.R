rm(list = ls())
source('./ecoFunctions.R')
source('./auxFunctions.R')
library(ggh4x)

# load data
files <- c(list.files('../data_sets', full.names = T),
           '../pyoverdine_data/training_set.csv', '../pyoverdine_data/test_set.csv')

data <- lapply(1:length(files),
               FUN = function(i) {
                 
                 df <- read.csv(files[i], header = T)
                 
                 if (i == 3) df$function. <- df$function. / 1e4
                 if (i %in% c(6, 7)) df <- cbind(df[, 1:8],
                                                 function. = rowMeans(df[, 9:ncol(df)]))
                 df <- aggregate(function. ~ .,
                                 data = df,
                                 FUN = mean)
                 
                 colnames(df)[ncol(df)] <- 'fun'
                 
                 return(df)
                 
               })
data[[6]] <- rbind(data[[6]], data[[7]])
data <- data[1:6]
data <- data[!grepl('butyrate', files[1:6])]

# evaluate predictions with/without residual inference
po <- do.call(rbind,
              lapply(data,
                     FUN = function(df) {
                       
                       po_res <- evaluatePredictions(df)
                       po_base <- evaluatePredictions_base(df)
                       
                       po_res$fun_predicted[po_res$fun_predicted < 0] <- 0
                       po_base$fun_predicted[po_base$fun_predicted < 0] <- 0
                       
                       R2_res <- 1 - sum((po_res$fun_true - po_res$fun_predicted)^2) / sum((po_res$fun_true - mean(po_res$fun_true))^2)
                       R2_base <- 1 - sum((po_base$fun_true - po_base$fun_predicted)^2) / sum((po_base$fun_true - mean(po_base$fun_true))^2)
                       
                       return(data.frame(R2_res = R2_res,
                                         R2_base = R2_base))
                       
                     }))

po$dataset <- basename(files[1:6])[!grepl('butyrate', files[1:6])]
po$dataset <- setNames(c('Bacterial starch hydrolysis',
                         'Phytoplankton biomass',
                         'Above-ground plant biomass',
                         'Bacterial xylose oxidation',
                         'Bacterial pyoverdine secretion'),
                       basename(files[1:6])[!grepl('butyrate', files[1:6])])[po$dataset]
po$dataset <- factor(po$dataset,
                     levels = c('Bacterial pyoverdine secretion',
                                'Above-ground plant biomass',
                                'Phytoplankton biomass',
                                'Bacterial xylose oxidation',
                                'Bacterial starch hydrolysis'))

po <- gather(po, method, R2, R2_res:R2_base)

# plot
mycolors <- c('black',
              '#d6d62d',
              '#66b666',
              '#cb96c3',
              '#d72027')

ggplot(po, aes(x = method, y = R2, color = dataset, group = dataset)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(labels = c('No residual\ninference',
                              'Residual inference')) +
  scale_y_continuous(name = expression(italic(R)^2~predictions~vs.~observations),
                     limits = c(0, 1)) +
  scale_color_manual(name = 'Data set',
                     values = mycolors) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 3,
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 14,
                                   angle = 90, hjust = 1),
        axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        axis.ticks.x = element_blank())

ggsave(filename ='../plots/res_vs_nores.pdf',
       device = 'pdf',
       dpi = 600,
       width = 150,
       height = 200,
       units = 'mm',
       limitsize = F)

