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
  
  species <- sort(unique(unlist(lapply(data[, 1],
                                       FUN = function(community) strsplit(community, split = ',')[[1]]))))
  
  data_out <- lapply(1:nrow(data),
                     FUN = function(i) matrix(as.numeric(species %in% strsplit(data[i, 1], split = ',')[[1]]), nrow = 1))
  data_out <- do.call(rbind, data_out)
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
if (run_tests) {
  test_that('Number of species',{
    expect_equal(nSpecies('1,2,3'), 3)
    expect_equal(nSpecies('lorem,ipsum'), 2)
    expect_equal(nSpecies(NA), 0)
  })
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
                              f_knockins <- f[names(f) %in% knockins]
                              f_knockins <- data.frame(community = names(f_knockins),
                                                       fun.knockin = as.numeric(f_knockins))
                              
                              f_backgrounds <- f[names(f) %in% backgrounds]
                              f_backgrounds <- data.frame(community = names(f_backgrounds),
                                                          fun.backgrund = as.numeric(f_backgrounds))
                              
                              f_knockins$background_community <- sapply(f_knockins$community,
                                                                        FUN = function(x) {
                                                                          xout <- strsplit(x, split = ',')[[1]]
                                                                          xout <- xout[xout != sp]
                                                                          xout <- paste(xout, collapse = ',')
                                                                          return(xout)
                                                                        })
                              
                              fun <- merge(f_backgrounds, f_knockins, by.x = 'community', by.y = 'background_community', all = T, suffixes = c('.bg', '.knockin'))
                              fun$df <- fun$fun.knockin - fun$fun.backgrund
                              fun <- fun[!is.na(fun$fun.backgrund) & !is.na(fun$fun.knockin), ] # keep only rows where the function of both the background and the knockin are known
                              
                              # build data frame
                              df <- data.frame(background = fun$community,
                                               knock_in = rep(sp, nrow(fun)),
                                               background_f = fun$fun.backgrund,
                                               d_f = fun$df)
                              
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

closestPaths <- function(target_comm, comm_list, species.removal = F) {
  
  # given a focal community (comm) and a list of other communities (comm_list), returns those in comm_list that are closest to comm (in terms of species addition)
  
  # fetch communities that are potential ancestors of the target
  isAncestor <- function(community_1, community_2) { # check if community_1 is an ancestor of community_2 (i.e. all species in community_1 are also in community_2)
    sapply(community_1,
           FUN = function(x) all(strsplit(x, split = ',')[[1]] %in% strsplit(community_2, split = ',')[[1]]))
  }
  if (!species.removal) comm_list <- comm_list[isAncestor(comm_list, target_comm)]
  
  if (length(comm_list) == 0) {
    
    warning(paste('For community\n', target_comm, '\nno ancestors in data. Returning <empty> community.', sep = ''))
    
    comm_list <- ''
    
  }
    
  all_comms <- data.frame(community = c(target_comm, comm_list),
                          fun = NA)
  all_comms <- string2matrix(all_comms)
  
  comm <- all_comms[1, 1:(ncol(all_comms) - 1), drop = F]
  comm_list <- all_comms[2:nrow(all_comms), 1:(ncol(all_comms) - 1), drop = F]
  
  dist <- sapply(1:nrow(comm_list),
                 FUN = function(i) sum(abs(comm_list[i, ] - comm)))
  
  which_min <- which(dist == min(dist))
  closest_comms <- comm_list[which_min, , drop = F]
    
  
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
               
               cumtraj_i <- sapply(cumtraj_i, # if there are repetitions of species, that means the species is being removed instead of added
                                   FUN = function(x) {
                                     
                                     xtab <- table(strsplit(x, split = ',')[[1]])
                                     removed_sp <- names(xtab)[xtab > 1]

                                     xout <- strsplit(x, split = ',')[[1]]
                                     xout <- paste(xout[!(xout %in% removed_sp)], collapse = ',')
                                     
                                     return(xout)
                                     
                                   })
               
               eps_i <- paste(cumtraj_i[-length(cumtraj_i)], traj_i, sep = '+')
               eps_i <- sapply(eps_i, # if knocked-in species is in background, it means it's a species removal
                               FUN = function(x) {
                                 
                                 knock_in <- strsplit(x, split = '\\+')[[1]][2]
                                 background <- strsplit(x, split = '\\+')[[1]][1]
                                 background <- strsplit(background, split = ',')[[1]]
                                 
                                 xout <- paste(background[background != knock_in], collapse = ',')
                                 xout <- paste(xout, knock_in, sep = '+')
                                 
                                 return(xout)
                                 
                               })
               eps_i <- as.character(eps_i)
               
               return(eps_i)
               
             })))
    
  }
  
  return(eps)
  
}

pathSteps <- function(target, source, single.traj = F) {

  # get steps in trajectory from source to target communities through addition/removal of species
  # is the closure condition is satisfied, a single trajectory is sufficient (speeds up calculations)
  
  path <- closestPaths(target, source)
  
  traj <- path[1, 4:ncol(path), drop = F]
  traj <- colnames(traj)[abs(traj[1, ]) == 1]
  
  # should a single trajectory be returned? (this is enough if the closer condition is satisfied)
  if (single.traj) {
    
    traj <- matrix(traj, nrow = 1)
    
  } else { # otherwise, return all possible trajectories (orders of species addition)
    
    if (length(traj) <= 4) {
      
      traj <- do.call(rbind, permn(traj))
      
    } else {
      
      # if number of steps is too large, getting all possible trajectories is computationally out of reach
      # in that case, we take only 25 arbitrarily chosen paths
      traj <- do.call(rbind,
                      lapply(1:25,
                             FUN = function(i) matrix(sample(traj, size = length(traj), replace = F), nrow = 1)))
      
    }
    
  }
  
  traj_steps <- lapply(1:nrow(traj),
                       FUN = function(i) {
                         
                         # trajectory
                         cumtraj <- orderName(sapply(1:length(traj[i, ]), FUN = function(j) paste(traj[i, 1:j], collapse = ',')))
                         cumtraj <- c(source,
                                      sapply(1:length(cumtraj), FUN = function(j) orderName(paste(c(strsplit(source, split = ',')[[1]],
                                                                                                    cumtraj[j]), collapse = ','))))
                         cumtraj <- sapply(cumtraj, # if there are repetitions of species, that means the species is being removed instead of added
                                           FUN = function(x) {
                                             
                                             xtab <- table(strsplit(x, split = ',')[[1]])
                                             removed_sp <- names(xtab)[xtab > 1]
                                             
                                             xout <- strsplit(x, split = ',')[[1]]
                                             xout <- paste(xout[!(xout %in% removed_sp)], collapse = ',')
                                               
                                             return(xout)
                                               
                                           })
                         cumtraj <- as.character(cumtraj)
                         
                         # knock-ins at each step
                         knockins <- traj[i, ]
                         
                         # backgrounds
                         backgrounds <- cumtraj[-length(cumtraj)]
                         backgrounds <- sapply(1:length(backgrounds),
                                               FUN = function(i) {
                                                 xout <- strsplit(backgrounds[i], split = ',')[[1]]
                                                 xout <- xout[xout != knockins[i]]
                                                 xout <- paste(xout, collapse = ',')
                                                 return(xout)
                                               })
                         backgrounds <- as.character(backgrounds)
                         
                         # signs
                         signs <- sapply(knockins, FUN = function(x) path[1, x])
                         signs <- as.numeric(signs)
                         
                         return(list(trajectory = cumtraj,
                                backgrounds = backgrounds,
                                knock_ins = knockins,
                                signs = signs))
                         
                       })
  
  if(single.traj) traj_steps <- traj_steps[[1]]
  
  return(traj_steps)

}

makeFEEs <- function(ge_data) {
  
  # make linear fits (input is a data.frame generated with makeGEdata)
  
  species <- sort(unique(ge_data$knock_in))
  fits <- do.call(rbind,
                  lapply(species,
                         FUN = function(sp){
                           data.frame(a = lm(data = ge_data[ge_data$knock_in == sp, ], formula = d_f~background_f)$coefficients[1],
                                      b = lm(data = ge_data[ge_data$knock_in == sp, ], formula = d_f~background_f)$coefficients[2])
                         }))
  rownames(fits) <- species
  
  return(fits)
  
}

inferAllResiduals <- function(ge_data) {

  # apply closure condition to the WHOLE landscape to infer values of residuals
  # only advisable for <10 species (otherwise pseudoinverse computation explodes)
  
  # get linear fits
  fits <- makeFEEs(ge_data)
  
  # fetch species names and build all combinatorial arrangements
  species <- sort(unique(ge_data$knock_in))
  communities <- orderName(unlist(sapply(1:length(species),
                                         FUN = function(i) sapply(combn(species,
                                                                        i,
                                                                        simplify = FALSE),
                                                                  FUN = paste, collapse = ','))))
  communities <- c('<empty>', communities) # attach the 'empty' community
  
  # list epsilons that are known from observations
  ge_data$background[ge_data$background == ''] <- '<empty>'
  known_eps <- setNames(ge_data$d_f- (fits[ge_data$knock_in, 'a'] + fits[ge_data$knock_in, 'b']*ge_data$background_f),
                        paste(ge_data$background, ge_data$knock_in, sep = '+'))
  
  # estimate standard deviations of epsilons
  sigma <- sapply(species,
                  FUN = function(sp) sd(known_eps[grepl(paste('\\+', sp, sep = ''), names(known_eps))]))
  
  # arrange epsilons in matrix form (rows: backgrounds, columns: knock-ins), this makes them easier to access later on
  eps_matrix <- matrix(NA, nrow = length(communities), ncol = length(species))
  colnames(eps_matrix) <- species
  rownames(eps_matrix) <- communities
  
  # some of the epsilons are known from the observations: add them to the matrix
  for (i in 1:length(known_eps)) {
    eps_matrix[strsplit(names(known_eps)[i], split = '\\+')[[1]][1],
               strsplit(names(known_eps)[i], split = '\\+')[[1]][2]] <- known_eps[i]
  }
  
  # elements of the matrix where the knock-in is already present in the background are set to NaN, other elements are either numeric (known epsilons) or NA (incognitas)
  for (sp in species) {
    eps_matrix[containsSpecies(sp, rownames(eps_matrix)), sp] <- NaN
  }
  
  # epsilons in vector form (makes it easier to operate with them)
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
  
  rownames(eps) <- gsub('<empty>', '', rownames(eps))
  communities[communities == '<empty>'] <- ''
  
  # build M matrix
  # i and j are two arbitrary species
  # B is an arbitrary background (not containing i or j), Bi and Bj are the same backgrounds containing also i and j respectively
  # a_i, b_i are the intercept and slope of the fit FEE for species i
  # e_i(B) is the epsilon for species i on background B
  # general form of the constraint:
  #   a_i*b_j + (1 + b_j)*e_i(B) + e_j(Bi) = a_j*b_i + (1 + b_i)*e_j(B) + e_i(Bj)
  constraints <- data.frame(B = character(0), i = character(0), j = character(0))
  
  # list all backgrounds with at least 2 species less than the most rich community (containing all species)
  bgs_valid <- communities[nSpecies(communities) <= length(species) - 2]
  
  # for each of those backgrounds, get all possible pairs of additional species (not in the background itself)
  sp_valid <- containsSpecies(species, bgs_valid)
  for (i in 1:length(bgs_valid)) {
    sp <- species[!sp_valid[, i]]
    pairs <- t(combn(sp, 2))
    constraints <- rbind(constraints,
                         data.frame(B = as.character(bgs_valid[i]),
                                    i = pairs[, 1],
                                    j = pairs[, 2]))
  }
  
  # some of the constraints are redundant; to identify them we need to check which edges of the 'community graph' are covered by each constraint
  isDescendant <- function(community_1, community_2) { # check if community_2 is a descendant of community_1 (i.e. all species in community_1 are also in community_2)
    sapply(community_2,
           FUN = function(x) all(strsplit(community_1, split = ',')[[1]] %in% strsplit(x, split = ',')[[1]]))
  }
  
  edges <- NULL
  nSpecies_comm <- nSpecies(communities)
  for (i in 1:length(communities)) {
    edges <- c(edges,
               paste(communities[i],
                     communities[nSpecies_comm == (1 + nSpecies_comm[i]) & isDescendant(communities[i], communities)],
                     sep = ' / '))
  }
  edges <- data.frame(edge = edges,
                      covered = FALSE)
  
  redundant_constraints <- NULL
  for (i in 1:nrow(constraints)) {
    
    # each constraint covers 4 communities
    c1 <- constraints[i, 1]
    c2 <- orderName(paste(c(strsplit(constraints[i, 1], split = ',')[[1]], constraints[i, 2]), collapse = ','))
    c3 <- orderName(paste(c(strsplit(constraints[i, 1], split = ',')[[1]], constraints[i, 3]), collapse = ','))
    c4 <- orderName(paste(c(strsplit(constraints[i, 1], split = ',')[[1]], constraints[i, 2], constraints[i, 3]), collapse = ','))
    
    # each constraint covers 4 edges
    edges_i <- c(paste(c1, c2, sep = ' / '),
                 paste(c1, c3, sep = ' / '),
                 paste(c2, c4, sep = ' / '),
                 paste(c3, c4, sep = ' / '))
    
    # which of those are already covered?
    covered_i <- edges$covered[edges$edge %in% edges_i]
    
    if (all(covered_i)) { # if they are all covered, tag the i-th constraint as redundant
      redundant_constraints <- c(redundant_constraints, i)
    } else { # if they are not all covered, keep the constraint and set them all to covered
      edges$covered[edges$edge %in% edges_i] <- TRUE
    }
    
  }
  
  # remove redundant constraints
  constraints <- constraints[-redundant_constraints, ]
  
  # formulate each constraint in matrix form
  M <- matrix(0, nrow = nrow(constraints), ncol = nrow(eps))
  C <- matrix(0, nrow = nrow(constraints), ncol = 1)
  rownames(M) <- paste('C.', 1:nrow(constraints), sep = '')
  rownames(C) <- rownames(M)
  colnames(M) <- rownames(eps)
  
  for (k in 1:nrow(M)) {
    
    i <- constraints$i[k]
    j <- constraints$j[k]
    B <- constraints$B[k]
    
    Bi <- orderName(paste(c(strsplit(B, split = ',')[[1]], i), collapse = ','))
    Bj <- orderName(paste(c(strsplit(B, split = ',')[[1]], j), collapse = ','))
    
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
    
    C[k] <- a_j*b_i - a_i*b_j
    
  }
  
  # check if any constraint is automatically satisfied by the known epsilons
  check <- M %*% eps - C
  check <- sapply(1:nrow(check), FUN = function(i) check[i] < 1e-5) # tolerance: 1e-5 (elements won't be exactly 0 even when the constraint is satisfied due to floating-point precision)
  
  stopifnot(all(check[!is.na(check)] == TRUE)) # all checks should either be TRUE (constraint satisfied) or NA (constraint not evaluated due to missing epsilons), if this is not the case something has failed
  
  # the dimension of M can be further reduced by moving the known epsilons to the right-hand side of the equation M %*% eps = C
  eps_1 <- eps
  eps_1[!is.na(eps)] <- 0
  eps_2 <- eps
  eps_2[is.na(eps_2)] <- 0 # eps = eps_1 + eps_2 ---> M %*% eps_1 = C - M %*% eps_2
  
  C_prime <- C - (M %*% eps_2) # this makes it so the equation to solve is M %*% eps_1 = C_prime
  remove_these <- sapply(1:nrow(M), FUN = function(i) all(as.numeric(M[i, ]) == rep(0, ncol(M)))) # to remove constraints that are automatically satisfied by the known epsilons
  C_prime <- C_prime[!remove_these, ]
  M <- M[!remove_these, 1:sum(is.na(eps))] # also removes the columns corresponding to known epsilons
  
  # scale elements of M by the sigmas
  M_sigma <- M * matrix(rep(sigma[gsub('.*\\+' ,'', colnames(M))], nrow(M)), nrow = nrow(M), byrow = TRUE)
  
  # now we just need to solve the system of linear equations M_sigma %*% eps_sigma = C_prime where the elements of eps_sigma are the incognitas
  # use the pseudoinverse for the solution that minimizes sum of squares of the elements of eps_sigma
  M_sigma_pseudoinv <- ginv(M_sigma)
  eps_sigma_0 <- M_sigma_pseudoinv %*% C_prime
  rownames(eps_sigma_0) <- colnames(M)
  
  # undo the sigma transformation to recover the epsilons
  eps_1 <- eps_sigma_0 * sigma[gsub('.*\\+' ,'', rownames(eps_sigma_0))]
  
  # add these values to the matrix of epsilons
  for (i in 1:nrow(eps_1)) {
    B <- gsub('\\+.*', '', rownames(eps)[i])
    if (B == '') B <- '<empty>'
    knockin <- gsub('.*\\+', '', rownames(eps)[i])
    eps_matrix[B, knockin] <- eps_1[i]
  }
  
  # return matrix of epsilons
  return(eps_matrix)

}

predictF_base <- function(target, data) {

  # predict the function of a community (target) starting from the function of an in-sample community
  # considers every unknown residual to be zero
  # first obtains the in-sample communities that are closer to the target community, then every possible path (i.e. order of addition/removal of species)
  # final prediction is the average across all paths
  
  ge_data <- makeGEdata(data)
  fits <- makeFEEs(ge_data)
  
  fun <- sapply(target,
                FUN = function(t) {
                  
                  paths <- closestPaths(t, data$community)
                  
                  # if there are too many possible paths from an in-sample community to the target community,
                  # we randomly choose 5 of them to avoid heavily increasing computation time
                  if (nrow(paths) > 5) paths <- paths[order(sample(1:nrow(paths), size = 5, replace = F)), ]
                  
                  fun <- lapply(1:nrow(paths),
                                FUN = function(i) {
                                  
                                  path <- paths[i, ]
                                  steps <- pathSteps(t, path$source[1], single.traj = F)
                                  
                                  fun_i <- sapply(steps,
                                                  FUN = function(st) {
                                                    
                                                    st$backgrounds[st$backgrounds == ''] <- '<empty>'
                                                    
                                                    # iterate
                                                    comm <- st$trajectory[1]
                                                    fun <- mean(data$fun[data$community == comm])
                                                    for (s in 1:length(st$knock_ins)) {
                                                      comm <- st$trajectory[s]
                                                      fun <- fun + st$signs[s]*(fits[st$knock_ins[s], 'a'] + fits[st$knock_ins[s], 'b']*fun)
                                                    }
                                                    
                                                    return(fun)
                                                    
                                                  })
                                  
                                  return(fun_i)
                                  
                                })
                  
                  return(mean(unlist(fun)))
                  
                })
  
  return(data.frame(community = names(fun),
                    fun = as.numeric(fun)))

}

predictF_fullClosure <- function(target, data, eps) {

  # predict the function of a community (target) starting from the function of an in-sample community
  # assuming that the closure condition is satisfied in the WHOLE landscape
  # therefore prediction should be independent from the path chosen (i.e. the source community and the order of species addition/removal)
  
  # format data just in case
  colnames(data) <- c('community', 'fun')
  data$community <- orderName(data$community)
  target <- orderName(target)
  
  # fit FEEs
  ge_data <- makeGEdata(data)
  fits <- makeFEEs(ge_data)
  
  # iterate over target (multiple target communities can be passed at once)
  predicted_f <- sapply(target,
                        FUN =  function(t) {
                          
                          # get closest path to target community (although any path should be equivalent if the closure condition is satisfied globally)
                          path <- closestPaths(t, data$community)
                          steps <- pathSteps(t, path$source[1], single.traj = T)
                          steps$backgrounds[steps$backgrounds == ''] <- '<empty>'
                          
                          # iterate
                          comm <- steps$trajectory[1]
                          fun <- mean(data$fun[data$community == comm])
                          for (s in 1:length(steps$knock_ins)) {
                            comm <- steps$trajectory[s]
                            fun <- fun + steps$signs[s]*(fits[steps$knock_ins[s], 'a'] + fits[steps$knock_ins[s], 'b']*fun + eps[steps$backgrounds[s], steps$knock_ins[s]])
                          }
                          
                          return(fun)
                          
                        })
  
  return(data.frame(community = names(predicted_f),
                    fun = as.numeric(predicted_f)))
  
}

if (F) {
  
  ### LOAD DATA FOR TESTING
  data <- read.table('../pyoverdine_data/pyo_rep3.txt', header = T)
  data <- rbind(data, data.frame(community = '', fun = 0))
  
  target <- sample(data$community, size = 30)
  if ('' %in% target) target <- target[target != '']
  obsF <- data[data$community %in% target, ]
  data <- data[!(data$community %in% target), ]
  
  predF <- predictF_base(target, data)
  
  po <- merge(obsF, predF, by = 'community', suffixes = c('_obs', '_pred'))
  plot(po$fun_pred, po$fun_obs)
  abline(a = 0, b = 1)

}
