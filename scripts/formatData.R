rm(list = ls())

# load auxiliary functions
source('./ecoFunctions.R')
library(scales)
library(gridExtra)

data <- lapply(list.files('../pyoverdine_data/', full.names = T), FUN = function(file) read.table(file, sep = '\t', header = T))
data <- merge(merge(data[[1]], data[[2]], all = T, by = 'community', suffixes = c('.rep1', '.rep2')),
              data[[3]], by = 'community')
colnames(data)[ncol(data)] <- 'fun.rep3'

data_training <- data[!is.na(data$fun.rep1), ]
data_test <- data[is.na(data$fun.rep1), c(1, 3, 4)]

data_training <- lapply(1:3,
                        FUN = function(i) string2matrix(data_training[, c(1, 1+i)]))

data_training <- merge(merge(data_training[[1]], data_training[[2]],
                             by = as.character(1:8), all = T, suffixes = c('.rep1', '.rep2')),
                       data_training[[3]],
                       by = as.character(1:8), all = T)
colnames(data_training) <- c(paste('sp_', as.character(1:8), sep = ''), paste('function.rep', 1:3, sep = ''))

data_test <- merge(string2matrix(data_test[, c(1,2)]), string2matrix(data_test[, c(1,3)]),
                   by = as.character(1:8), all = T, suffixes = c('.rep1', '.rep2'))
colnames(data_test) <- c(paste('sp_', as.character(1:8), sep = ''), paste('function.rep', 1:2, sep = ''))

write.csv(data_training, file = '../pyoverdine_data/training_set.csv', row.names = F, quote = F)
write.csv(data_test, file = '../pyoverdine_data/test_set.csv', row.names = F, quote = F)