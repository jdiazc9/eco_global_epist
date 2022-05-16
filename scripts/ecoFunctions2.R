rm(list = ls())

### LIBRARIES

require(testthat)
require(MASS)
require(ggplot2)
require(gtools)
require(scales)
require(gridExtra)
require(tidyverse)
require(combinat)

run_tests <- FALSE # set to TRUE to run tests in every execution

### FUNCTIONS & TESTS

orderName <- function(community) {
  
  # order community names with the form 'sp2,sp2,sp3,...' so that species consistently
  # appear in alphabetical order
  
  as.character(sapply(community,
                      FUN = function(s) paste(sort(strsplit(s, split = ',')[[1]]),
                                              collapse = ',')))
}

string2matrix <- function(data) {
  
  # takes a data frame with 2 columns, the first one containing community names
  # (species present separated by commas with no spaces, e.g. 'sp1,sp2,sp3,...') and
  # the second one community functions, and converts it to binary matrix format
  
  species <- unique(unlist(sapply(data[, 1],
                                  FUN = function(community) strsplit(community, split = ',')[[1]])))
  
  data_out <- t(sapply(1:nrow(data),
                       FUN = function(i) as.numeric(species %in% strsplit(data[i, 1], split = ',')[[1]])))
  colnames(data_out) <- species
  
  data_out <- cbind(data_out, fun = data[, 2])
  
  return(data_out)
  
}

matrix2string <- function(data) {
  
  # inverse operation with respect to string2matrix
  
  species <- colnames(data)[-ncol(data)]
  
  communities <- sapply(1:nrow(data),
                        FUN = function(i) paste(species[data[i, -ncol(data)] == 1], collapse = ','))
  communities <- orderName(communities)
  
  return(data.frame(community = communities,
                    fun = data[, ncol(data)]))
  
}

if (run_tests) {
  test_that('Name ordering',{
    expect_equal(orderName('1,3,2'), '1,2,3')
    expect_equal(orderName(c('B,C,A,F,D', 'y,x,z')), c('A,B,C,D,F', 'x,y,z'))
  })
}

nSpecies <- function(community) {
  
  # returns the number of species in a community (given in 'sp1,sp2,sp3,...' format)
  # returns 0 if input is NA
  
  as.numeric(sapply(community,
                    FUN = function(comm) {
                      if (is.na(comm)) return(0)
                      else return(length(strsplit(comm, split = ',')[[1]]))
                    }))
  
}

containsSpecies <- function(species, community) {
  
  # returns TRUE if species is present in community, FALSE otherwise
  
  return(sapply(community,
                FUN = function(comm) species %in% strsplit(comm, split = ',')[[1]]))
  
}

makeGEdata <- function(data) {
  
  # takes a mapping between community structures and functions and returns a data frame of species functional effects
  # input data must be a data frame where he first column corresponds to genotype names and the second column corresponds to genotype fitness
  # genotype names should be such that present mutations appear separated by commas with no spaces, e.g. 'sp1,sp2,sp3,...'
  # use matrix2string() function to convert data to this format if needed
  
  # name columns of data
  colnames(data) <- c('community', 'f')
  
  # order community names alphabetically
  data[, 1] <- orderName(data[, 1])
  
  # if there are multiple instances of a same genotype, take the average and print a warning
  if(!all(table(data[, 1]) == 1)) {
    warning(paste('Multiple instances of a same combination in input data set:\n',
                  paste(names(table(data[, 1]))[table(data[, 1]) > 1], collapse = '\n'),
                  '\nAveraging F to proceed.', sep = ''))
    data <- aggregate(formula = f ~ community,
                      data = data,
                      FUN = mean)
  }
  
  # extract species names
  species <- sort(unique(unlist(strsplit(data[, 1], split = ','))))
  n_species <- length(species)
  
  # functions as named array
  f <- setNames(data[, 2], data[, 1])
  
  # output: data frame where the effect of each species on a background community is isolated
  ge_data <- do.call(rbind,
                     lapply(species,
                            FUN = function(sp) {
                              
                              # fetch all communities in sample that contain species sp
                              knockins <- data$community[containsSpecies(sp, data$community)]
                              
                              # backgrounds corresponding to those communities
                              backgrounds <- as.character(sapply(knockins,
                                                                 FUN = function(x) {
                                                                   x <- strsplit(x, split = ',')[[1]]
                                                                   x <- x[x != sp]
                                                                   x <- paste(x, collapse = ',')
                                                                   return(x)
                                                                 }))
                              
                              # functions of backrounds and backgrounds + knock-ins
                              f_knockins <- as.numeric(f[knockins])
                              f_backgrounds <- as.numeric(f[backgrounds])
                              
                              # build data frame
                              df <- data.frame(background = backgrounds,
                                               knock_in = sp,
                                               background_f = f_backgrounds,
                                               d_f = f_knockins - f_backgrounds)
                              
                              # keep only rows where the function of both the background and the knockin are known
                              df <- df[!is.na(df$background_f) & !is.na(df$d_f), ]
                              
                              rownames(df) <- NULL
                              
                              return(df)
                              
                            }))
  
  return(ge_data)
  
}

plotFEEs <- function(ge_data) {
  
  # make dF-vs-F plots
  
  p <-
    ggplot(ge_data, aes(x = background_f, y = d_f)) +
      geom_abline(slope = 0,
                  intercept = 0,
                  color = '#d1d3d4') +
      geom_point(shape = 1,
                 cex = 2) +
      geom_smooth(method = 'lm',
                  formula = y~x,
                  color = 'firebrick1',
                  se = F,
                  fullrange = T) +
      scale_x_continuous(breaks = pretty_breaks(n = 3),
                         name = 'Function of ecological background [a.u.]') +
      scale_y_continuous(breaks = pretty_breaks(n = 3),
                         name = 'dF [a.u.]') +
      facet_wrap(~knock_in) +
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
  
  print(p)
  
  return(p)
  
}

closestPaths <- function(target_comm, comm_list) {
  
  # given a focal community (comm) and a list of other communities (comm_list), returns those in comm_list that are closest to comm (in terms of species addition/removal)
  
  all_comms <- data.frame(community = c(target_comm, comm_list),
                          fun = NA)
  all_comms <- string2matrix(all_comms)
  
  comm <- all_comms[1, 1:(ncol(all_comms) - 1)]
  comm_list <- all_comms[2:nrow(all_comms), 1:(ncol(all_comms) - 1)]
  
  dist <- sapply(1:nrow(comm_list),
                 FUN = function(i) sum(abs(comm_list[i, ] - comm)))
  
  which_min <- which(dist == min(dist))
  closest_comms <- comm_list[which_min, ]
  
  if (length(which_min) == 1) {
    closest_comms <- matrix(closest_comms, nrow = 1)
    colnames(closest_comms) <- names(comm)
  }
  
  # paths to comm
  paths <- lapply(1:nrow(closest_comms),
                  FUN = function(i) comm - closest_comms[i, ])
  paths <- do.call(rbind, paths)
  
  sources <- matrix2string(cbind(as.data.frame(closest_comms), fun = NA))$community
  
  # return info
  return(cbind(source = sources, target = target_comm, dist = dist[which_min], as.data.frame(paths)))
  
}

fetchResiduals <- function(paths) {
  
  # fetch which residuals (epsilons) take part in a given set of paths (from a source to a target community)
  # accounts for every possible path order (e.g. i->j->k, j->i->k, k->i->j, ...)
  # variable naming: eps_i(s) is named s+i (e.g. 'i,j+k', 'i,k+j', 'i,j,k+l', ...)
  
  eps <- NULL
  
  for (p in 1:nrow(paths)) {
    
    traj <- paths[p, 4:ncol(paths)]
    traj <- colnames(traj)[abs(traj[1, ]) == 1]
    traj <- do.call(rbind, permn(traj))
    
    eps <- c(eps,
             unlist(lapply(1:nrow(traj), FUN = function(t) {
               
               traj_i <- traj[t, ]
               cumtraj_i <- orderName(sapply(1:length(traj_i), FUN = function(i) paste(traj_i[1:i], collapse = ',')))
               cumtraj_i <- c(paths$source[p],
                              orderName(paste(paths$source[p], cumtraj_i, sep = ',')))
               eps_i <- paste(cumtraj_i[-length(cumtraj_i)], traj_i, sep = '+')
               
               return(eps_i)
               
             })))
    
  }
  
  return(eps)
  
}




### LOAD DATA FOR TESTING
data <- read.table('../pyoverdine_data/pyo_rep3.txt', header = T)
data <- rbind(data, data.frame(community = '', fun = 0))

ge_data <- makeGEdata(data)


set.seed(0)
target_comm <- data$community[100]
data <- data[-100, ]
data <- data[sample(1:nrow(data), size = 10), ]

paths <- closestPaths(target_comm, data$community)
paths2 <- paths
paths2$source <- '1,8'
paths2[, 4:11] <- c(0, 1, 1, 0, 0, 0, 0, 0)
paths <- rbind(paths, paths2)
rm(paths2)

fetchResiduals(paths)

























if(F) {




















nMut <- function(genotypes) {
  
  as.numeric(sapply(genotypes,
                    FUN = function(genotype) {
                      if (is.na(genotype)) return(0)
                      else return(length(strsplit(genotype, split = ',')[[1]]))
                    }))
}

if (run_tests) {
  test_that('Number of mutations',{
    expect_equal(nMut('1,2,3'), 3)
    expect_equal(nMut('lorem,ipsum'), 2)
    expect_equal(nMut(NA), 0)
  })
}












# format data: from an (incomplete) genotype-to-fitness map to 'species fitness effects' data
makeGEdata <- function(data, exclude.single.mut = FALSE, baseline.fun = 0) {
  
  # data must be a data frame where he first column corresponds to genotype names and the second column corresponds to genotype fitness
  # genotype names should be such that present mutations appear separated by commas with no spaces, e.g. 'i,j,k,...'
  
  # name columns of data
  colnames(data) <- c('genot', 'f')
  
  # order community names alphabetically
  data[, 1] <- orderNames(data[, 1])
  
  # if there are multiple instances of a same genotype, take the average and print a warning
  if(!all(table(data[, 1]) == 1)) {
    warning(paste('Multiple instances of a same combination in input data set:\n',
                  paste(names(table(data[, 1]))[table(data[, 1]) > 1], collapse = '\n'),
                  '\nAveraging F to proceed.', sep = ''))
    data <- aggregate(formula = f ~ genot,
                      data = data,
                      FUN = mean)
  }
  
  # extract mutation names
  muts <- sort(unique(unlist(strsplit(data[, 1], split = ','))))
  n_muts <- length(muts)
  
  # extract single mutant data
  single_mut <- data[nMut(data$genot) == 1, ]
  
  # fitness as named array
  f <- setNames(data[, 2], data[, 1])
  
  # output: data frame where the effect of each species on a background community is isolated
  ge_data <- do.call(rbind,
                     lapply(muts,
                            FUN = function(mut_i) {
                              
                              # all communities that contain mut_i
                              k_i <- data$genot[grepl(mut_i, data$genot)]
                              
                              # corresponding backgrounds
                              bg_i <- as.character(sapply(k_i,
                                                          FUN = function(x) {
                                                            x <- strsplit(x, split = ',')[[1]]
                                                            x <- x[x != mut_i]
                                                            x <- paste(x, collapse = ',')
                                                            return(x)
                                                          }))
                              
                              # functions of backrounds and backgrounds + knock-ins
                              f_k_i <- as.numeric(f[k_i])
                              f_bg_i <- as.numeric(f[bg_i])
                              
                              # return as data frame
                              ge_data_i <- data.frame(background = bg_i,
                                                      knock_in = mut_i,
                                                      background_f = f_bg_i,
                                                      d_f = f_k_i - f_bg_i)
                              ge_data_i <- ge_data_i[!is.na(ge_data_i$background_f) & !is.na(ge_data_i$d_f), ]
                              rownames(ge_data_i) <- NULL
                              
                              return(ge_data_i)
                              
                            }))
  
  # attach monoculture data if required, background community names set to NA and functions set to baseline.fun (0 by default)
  # (could be set to 1 if e.g. the 'no background' instance corresponds to a wild-type genotype)
  if(nrow(single_mut) > 0 & exclude.single.mut == FALSE) {
    ge_data <- rbind(ge_data,
                     data.frame(background = NA,
                                knock_in = single_mut$genot,
                                background_f = baseline.fun,
                                d_f = single_mut$f))
  }
    
  return(ge_data)
  
}

# plot global epistasis patterns (expects data in GE format)
plotGE <- function(ge_data, title = '') {
  
  # If raw.data = TRUE, this function expects the data in the global input format,
  # i.e. a two-column data frame where the first column corresponds to genotype names
  # (with mutations separated by commas and no spaces: A,B,C,...)
  # and the second column contains fitness values for each genotype.
  # In this case, data are first formatted through the makeGEdata function.
  # If raw.data = FALSE, the function expects a data frame that has already been
  # formatted through the makeGEdata function.
  
  geplot <- 
    ggplot(ge_data, aes(x = background_f, y = d_f, color = knock_in)) +
    geom_abline(slope = 0,
                intercept = 0,
                linetype = 'dashed') +
    geom_point() +
    geom_smooth(method = 'lm',
                formula = y~x,
                se = FALSE) +
    scale_x_continuous(name = 'F (background)') +
    scale_y_continuous(name = 'dF') +
    facet_wrap(~knock_in) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          legend.position = 'none') +
    ggtitle(title)
  
  print(geplot)
  return(geplot)
  
}

# make linear fits assuming no branching (expects data in GE format)
makeGEfits <- function(ge_data) {
  
  muts <- unique(ge_data$knock_in)
  fits <- do.call(rbind,
                  lapply(muts,
                         FUN = function(mut_i){
                           data.frame(a = lm(data = ge_data[ge_data$knock_in == mut_i, ], formula = d_f~background_f)$coefficients[1],
                                      b = lm(data = ge_data[ge_data$knock_in == mut_i, ], formula = d_f~background_f)$coefficients[2])
                         }))
  rownames(fits) <- muts
  
  return(fits)
  
}

# infer the values of the unknown epsilons (expects data in GE format)
inferEps <- function(ge_data) {
  
  # split single mutants from rest of the data
  single_mut <- ge_data[nMut(ge_data$background) == 0, ]
  ge_data <- ge_data[nMut(ge_data$background) > 0, ]
  
  # get linear fits
  fits <- makeGEfits(ge_data)
  
  # fetch mutation names and all combinatorial landscape
  muts <- unique(ge_data$knock_in)
  genots <- orderNames(unlist(sapply(1:length(muts),
                                     FUN = function(i) sapply(combn(muts,
                                                                    i,
                                                                    simplify = FALSE),
                                                              FUN = paste, collapse = ','))))
  
  # list epsilons that are known from observations
  known_eps <- setNames(ge_data$d_f- (fits[ge_data$knock_in, 'a'] + fits[ge_data$knock_in, 'b']*ge_data$background_f),
                        paste(ge_data$background, ge_data$knock_in, sep = '+'))
  
  # estimate standard deviations of epsilons
  sigma <- sapply(muts,
                  FUN = function(mut_i) sd(known_eps[grepl(paste('+', mut_i, sep = ''), names(known_eps))]))
  
  # arrange epsilons in matrix form (rows: backgrounds, columns: knock-ins), this makes them easier to access later on
  eps_matrix <- matrix(NA, nrow = length(genots), ncol = length(muts))
  colnames(eps_matrix) <- muts
  rownames(eps_matrix) <- genots
  for (i in 1:length(known_eps)) {
    eps_matrix[strsplit(names(known_eps), split = '\\+')[[i]][1],
               strsplit(names(known_eps), split = '\\+')[[i]][2]] <- known_eps[i]
  }
  for (mut_i in muts) {
    eps_matrix[grepl(mut_i, genots), mut_i] <- NaN # elements of the matrix where the knock-in is already present in the background are set to NaN, other elements are either numeric (known epsilons) or NA (incognitas)
  }
  
  # epsilons in vector form
  eps <- as.data.frame(matrix(NA, nrow = length(eps_matrix), ncol = 2))
  for (i in 1:nrow(eps_matrix)) {
    for (j in 1:ncol(eps_matrix)) {
      eps[ncol(eps_matrix)*(i-1) + j, 1] <- paste(rownames(eps_matrix)[i], colnames(eps_matrix)[j], sep = '+')
      eps[ncol(eps_matrix)*(i-1) + j, 2] <- eps_matrix[i, j]
    }
  }
  eps <- eps[!is.nan(eps[, 2]), ]
  eps <- eps[order(eps[, 2], na.last = FALSE), ]
  rownames(eps) <- eps[, 1]
  eps <- matrix(eps[, 2], dimnames = list(eps[, 1], 'eps'))
  
  # build M matrix
  # i and j are two arbitrary mutations
  # B is an arbitrary background (not containing i or j), Bi and Bj are the same backgrounds containing also i and j respectively
  # a_i, b_i are the intercept and slope of the fit for species i
  # e_i(B) is the epsilon for species i on background B
  # general form of the constraint:
  #   a_i*b_j + (1 + b_j)*e_i(B) + e_j(Bi) = a_j*b_i + (1 + b_i)*e_j(B) + e_i(Bj)
  constraints <- data.frame(B = character(0), i = character(0), j = character(0))
  
  # list all backgrounds with at least 2 mutations less than the 'all-mutated' genotype
  bgs_valid <- genots[nMut(genots) < length(muts) - 1]
  
  # for each of those backgrounds, get all possible pairs of additional mutations (not in the background itself)
  for (i in 1:length(bgs_valid)) {
    mut_i <- muts[sapply(muts, FUN = function(x) !(x %in% strsplit(bgs_valid[i], split = ',')[[1]]))]
    pairs_i <- t(combn(mut_i, 2))
    constraints <- rbind(constraints,
                         data.frame(B = as.character(bgs_valid[i]),
                                    i = pairs_i[, 1],
                                    j = pairs_i[, 2]))
  }
  
  # some of the constraints are redundant; to identify them we need to check which edges of the 'fitness graph' are covered by each constraint
  isDescendant <- function(genot_1, genot_2) { # check if genot_2 is a descendant of genot_1
    sapply(genot_2,
           FUN = function(x) all(strsplit(genot_1, split = ',')[[1]] %in% strsplit(x, split = ',')[[1]]))
  }
  
  edges <- NULL
  nMut_genots <- nMut(genots)
  for (i in 1:length(genots)) {
    edges <- c(edges,
               paste(genots[i],
                     genots[nMut_genots == (1 + nMut_genots[i]) & isDescendant(genots[i], genots)],
                     sep = ' / '))
  }
  
  edges <- data.frame(edge = edges,
                      covered = FALSE)
  
  # check which edges are covered by each constraint; if a constraint covers only edges that are already covered, tag it as redundant
  redundant_constraints <- NULL
  for (i in 1:nrow(constraints)) {
    
    # each constraint covers 4 edges
    edges_i <- c(paste(constraints[i, 1],
                       orderNames(paste(constraints[i, 1], constraints[i, 2], sep = ',')),
                       sep = ' / '),
                 paste(constraints[i, 1],
                       orderNames(paste(constraints[i, 1], constraints[i, 3], sep = ',')),
                       sep = ' / '),
                 paste(orderNames(paste(constraints[i, 1], constraints[i, 2], sep = ',')),
                       orderNames(paste(constraints[i, 1], constraints[i, 2], constraints[i, 3], sep = ',')),
                       sep = ' / '),
                 paste(orderNames(paste(constraints[i, 1], constraints[i, 3], sep = ',')),
                       orderNames(paste(constraints[i, 1], constraints[i, 2], constraints[i, 3], sep = ',')),
                       sep = ' / '))
    
    # which of those are already covered?
    covered_i <- edges$covered[edges$edge %in% edges_i]
    
    # if they are all covered, tag the i-th constraint as redundant
    if (all(covered_i)) {
      redundant_constraints <- c(redundant_constraints, i)
    } else { # if they are not all covered, keep the constraint and set them all to covered
      edges$covered[edges$edge %in% edges_i] <- TRUE
    }
    
  }
  
  # get rid of redundant constraints
  constraints <- constraints[-redundant_constraints, ]
  
  # formulate each constraint in matrix form
  M <- matrix(0, nrow = nrow(constraints), ncol = nrow(eps))
  c <- matrix(0, nrow = nrow(constraints), ncol = 1)
  rownames(M) <- paste('C.', 1:nrow(constraints), sep = '')
  rownames(c) <- rownames(M)
  colnames(M) <- rownames(eps)
  
  for (k in 1:nrow(M)) {
    
    i <- constraints$i[k]
    j <- constraints$j[k]
    B <- constraints$B[k]
    
    Bi <- orderNames(paste(B, i, sep = ','))
    Bj <- orderNames(paste(B, j, sep = ','))
    
    a_i <- fits[i, 'a']
    b_i <- fits[i, 'b']
    a_j <- fits[j, 'a']
    b_j <- fits[j, 'b']
    
    e_i_B <- paste(B, i, sep = '+')
    e_j_B <- paste(B, j, sep = '+')
    e_j_Bi <- paste(Bi, j, sep = '+')
    e_i_Bj <- paste(Bj, i, sep = '+')
    
    M[k, e_i_B] <- 1 + b_j
    M[k, e_j_B] <- -(1 + b_i)
    M[k, e_j_Bi] <- 1
    M[k, e_i_Bj] <- -1
    
    c[k] <- a_j*b_i - a_i*b_j
    
  }
  
  # check if any constraint is automatically satisfied by the known epsilons
  check <- M %*% eps - c
  check <- sapply(1:nrow(check), FUN = function(i) check[i] < 1e-5) # tolerance: 1e-5 (elements won't be exactly 0 even when the constraint is satisfied due to floating-point precision)
  
  stopifnot(all(check[!is.na(check)] == TRUE)) # all checks should either be TRUE (constraint satisfied) or NA (constraint not evaluated due to missing epsilons), if this is not the case return a warning
  
  # leave only unknown epsilons as incognitas
  eps_1 <- eps
  eps_1[!is.na(eps)] <- 0
  eps_2 <- eps
  eps_2[is.na(eps_2)] <- 0 # eps = eps_1 + eps_2 ---> M %*% eps_1 = c - M %*% eps_2
  
  c_prime <- c - (M %*% eps_2) # this makes it so the equation to solve is M %*% eps_1 = c_prime
  remove_these <- sapply(1:nrow(M), FUN = function(i) all(as.numeric(M[i, ]) == rep(0, ncol(M)))) # to remove constraints that are automatically satisfied by the known epsilons
  c_prime <- c_prime[!remove_these, ]
  M <- M[!remove_these, 1:sum(is.na(eps))] # alos removes the columns corresponding to known epsilons
  
  # rewrite problem in terms of epsilon_sigma and M_sigma
  M_sigma <- M * matrix(rep(sigma[gsub('.*\\+' ,'', colnames(M))], nrow(M)), nrow = nrow(M), byrow = TRUE)
  
  # now we just need to solve the system of linear equations sigma_M %*% eps_sigma = c_prime where the elements of eps_sigma are the incognitas
  # FIXME: what is the dimension of M? after removing redundant constraints? how do we remove redundant constraints in practice? otherwise computation time might explode for large systems
  # use the pseudoinverse for the solution that minimizes sum of squares of the elements of eps_sigma
  M_sigma_pseudoinv <- ginv(M_sigma)
  eps_sigma_0 <- M_sigma_pseudoinv %*% c_prime
  rownames(eps_sigma_0) <- colnames(M)
  
  # undo the sigma transformation to recover the epsilons
  eps <- eps_sigma_0 * sigma[gsub('.*\\+' ,'', rownames(eps_sigma_0))]
  
  # add these values to the matrix of epsilons
  for (i in 1:nrow(eps)) {
    eps_matrix[gsub('\\+.*', '', rownames(eps)[i]),
               gsub('.*\\+', '', rownames(eps)[i])] <- eps[i]
  }
  
  # return matrix of epsilons
  return(eps_matrix)
  
}

predictF <- function(genot, data, coeff, eps) {
  
  # genot: the name(s) of the genotype(s) to predict, e.g. c('i,j,k', 'l,m')
  # single.mut: a data.frame containing the names (first column) and fitness (second column) of the single mutants
  # coeff: a data.frame of coefficients for the linear fits (obtained from makeGEfits)
  # eps: a matrix of epsilons (inferred from inferEps)
  
  data[, 1] <- orderNames(data[, 1])
  genot <- orderNames(genot)
  
  # aux. function: get the mutational path between two genotypes: what mutations are needed to reach genot_2 from genot_1?
  mutPath <- function(genot_1, genot_2) {
    
    genot_1 <- strsplit(genot_1, split = ',')[[1]]
    genot_2 <- strsplit(genot_2, split = ',')[[1]]
    
    lose_these <- genot_1[!(genot_1 %in% genot_2)]
    gain_these <- genot_2[!(genot_2 %in% genot_1)]
    
    return(list(path_length = length(gain_these) + length(lose_these),
                lose_these = lose_these,
                gain_these = gain_these))
    
  }
  
  predicted_f <- sapply(genot,
                        FUN = function(genot_i) {
                          
                          if(genot_i %in% data[, 1]) {
                            
                            return(data[data[, 1] == genot_i, 2])
                            
                          } else {
                          
                            # identify closest genotypes in dataset
                            path_lengths <- sapply(data[, 1], FUN = function(x) mutPath(x, genot_i)$path_length)
                            closest_genots <- names(path_lengths)[path_lengths == min(path_lengths)]
                            
                            # how many unknown epsilons are needed to reach genot_i from each of the closest_genots?
                            best_paths <- vector(mode = 'list', length = length(closest_genots))
                            names(best_paths) <- closest_genots
                            
                            for (i in 1:length(closest_genots)) {
                              
                              path <- mutPath(closest_genots[i], genot_i)
                              all_paths <- permutations(path$path_length, path$path_length,
                                                        v = c(path$lose_these, path$gain_these))
                              
                              n_eps <- rep(0, nrow(all_paths))
                              
                              for (j in 1:nrow(all_paths)) {
                                
                                this_step <- closest_genots[i]
                                
                                for (step in 1:path$path_length) {
                                  
                                  next_step <- this_step
                                  if (all_paths[j, step] %in% path$lose_these) {
                                    next_step <- strsplit(next_step, split = ',')[[1]]
                                    next_step <- next_step[next_step != all_paths[j, step]]
                                    next_step <- paste(next_step, collapse = ',')
                                  } else if (all_paths[j, step] %in% path$gain_these) {
                                    next_step <- orderNames(paste(next_step, all_paths[j, step], collapse = ','))
                                  }
                                  
                                  is_eps_known <- this_step %in% data[, 1] & next_step %in% data[, 1]
                                  n_eps[j] <- n_eps[j] + !is_eps_known
                                  
                                  this_step <- next_step
                                  
                                }
                                
                              }
                              
                              # get order with least number of unknown epsilons required (if there is a tie, take the first one)
                              path$n_eps <- min(n_eps)
                              path$order <- all_paths[which(n_eps == min(n_eps))[1], ]
                              
                              best_paths[[i]] <- path
                              
                            }
                            
                            ### FIXME: can we generalize this?
                            # discard paths that involve passing through the wild-type genotype (no mutations)
                            # currently we are not considering the epsilons of the form eps_i_B where B is the 'empty' background
                            discard_these <- sapply(1:length(best_paths),
                                                    FUN = function(p) {
                                                      
                                                      path_i <- names(best_paths)[p]
                                                      
                                                      for (i in 1:best_paths[[p]]$path_length) {
                                                        
                                                        path_next <- path_i[length(path_i)]
                                                        
                                                        if (best_paths[[p]]$order[i] %in% best_paths[[p]]$gain_these) {
                                                          path_next <- paste(path_next, best_paths[[p]]$order[i], sep = ',')
                                                        } else if (best_paths[[p]]$order[i] %in% best_paths[[p]]$lose_these) {
                                                          path_next <- strsplit(path_next, split = ',')[[1]]
                                                          path_next <- path_next[!(path_next == best_paths[[p]]$order[i])]
                                                          path_next <- paste(path_next, collapse = ',')
                                                        }
                                                        
                                                        path_i <- c(path_i, path_next)
                                                        
                                                      }
                                                      
                                                      return(any(path_i == ''))
                                                      
                                                    })
                            
                            best_paths[which(discard_these)] <- NULL
                            closest_genots <- closest_genots[!discard_these]
                            
                            # which path involves the least amount of unknown epsilons? (if there is a tie, take the first one)
                            which_best_path <- which.min(sapply(best_paths, FUN = function(x) x$n_eps))[1]
                            best_path <- best_paths[[which_best_path]]
                            best_path$source <- closest_genots[which_best_path]
                            
                            # follow the best path to predict function
                            fun <- data[data[, 1] == best_path$source, 2]
                            B <- best_path$source
                            
                            for (step in 1:best_path$path_length) {
                              
                              i <- best_path$order[step]
                              
                              if (best_path$order[step] %in% best_path$gain_these) {
                                
                                d_fun <- coeff[i, 'a'] + coeff[i, 'b']*fun + eps[B, i]
                                fun <- fun + d_fun
                                
                                B <- orderNames(paste(B, i, sep = ','))
                                
                              } else if (best_path$order[step] %in% best_path$lose_these) {
                                
                                B <- strsplit(B, split = ',')[[1]]
                                B <- paste(B[B != i], collapse = ',')
                                
                                d_fun <- coeff[i, 'a'] + coeff[i, 'b']*fun + eps[B, i]
                                fun <- fun - d_fun
                                
                              }
                              
                            }
                            
                            return(fun)
                          
                          }
                          
                        })
  
  return(predicted_f)
  
}




### ----------------------------------------------------------------------------

if (F) {


### RUN PIPELINE

# load data
data_full <- read.table('../jack/data/ge_pyoverdine.txt', header = T, sep = '\t')
data_full[, 1] <- orderNames(data_full[, 1])

# aggregate replicate communities if necessary (makeGEdata does this automatically but it's better to do it before hand for subsampling)
data_full <- aggregate(formula = fun ~ community,
                       data = data_full,
                       FUN = mean)

# GE formatted data
ge_data_full <- makeGEdata(data_full)

# fetch mutation names and build full combinatorial space
muts <- unique(ge_data_full$knock_in)
genots <- orderNames(unlist(sapply(1:length(muts),
                                   FUN = function(i) sapply(combn(muts,
                                                                  i,
                                                                  simplify = FALSE),
                                                            FUN = paste, collapse = ','))))

# predict entire landscape
fits <- makeGEfits(ge_data_full)
eps <- inferEps(ge_data_full)
predicted_f <- predictF(genots,
                        data_full,
                        fits,
                        eps)

# partial observations
observed_f <- setNames(data_full$fun, data_full$community)
observed_f <- observed_f[genots]
po <- data.frame(comunity = genots,
                 predicted_f = predicted_f,
                 observed_f = observed_f)
po <- po[order(po$predicted_f, decreasing = TRUE), ]
po$predicted_rank <- 1:nrow(po)

# subsample data (make sure that monoculture data is left in the sample)
leave_out_of_sample <- sample(data_full$community, size = round(0.2*nrow(data_full)), replace = FALSE)
leave_out_of_sample <- leave_out_of_sample[nMut(leave_out_of_sample) > 1]

data <- data_full[!(data_full$community %in% leave_out_of_sample), ] # subsample data

# pipeline
geplot <- plotGE(data)
fits <- makeGEfits(data)
eps <- inferEps(data)
predicted_f <- predictF(leave_out_of_sample,
                        data,
                        fits,
                        eps)
observed_f <- data_full[data_full$community %in% leave_out_of_sample, ]

pred_obs_f <- merge(observed_f,
                    data.frame(community = leave_out_of_sample,
                               fun = predicted_f),
                    by = 'community',
                    suffixes = c('_obs', '_pred'))

ggplot(pred_obs_f, aes(x = fun_pred, y = fun_obs)) +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed') +
  geom_point() +
  geom_smooth(method = 'lm',
              se = FALSE,
              formula = y~x,
              color = 'black') +
  scale_x_continuous(name = 'Predicted F') +
  scale_y_continuous(name = 'Observed F') +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

print(cor(pred_obs_f$fun_obs, pred_obs_f$fun_pred)^2)

}



}
