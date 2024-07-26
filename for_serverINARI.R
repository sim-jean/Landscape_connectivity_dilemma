### Generate small scale landscapes on INARI ####

rm(list = ls())

library(Matrix)
library(parallel)
if(!("digest" %in% installed.packages())){
  install.packages("digest")
}
library(digest)

## 
print("Environment is set up ok")

file_path_ = "/home/jean/connectivity_dilemma/"

numCores = 10

values = c(0,1,2)


for(size in 3:4){
  
  print(paste("Size", size, 'is starting'))
  equi_landscape = function(land, include_ = F){
    # 3 transformations clockwise rotation
    m2  =  t(apply(land, 2, rev))
    m3  =  t(apply(m2, 2, rev))
    m4  =  t(apply(m3, 2, rev))
    
    # Horizontal symmetry
    hor_sym  =  land[nrow(land):1, ]
    
    # 3 clockwise transformations of horizontally symmetrized matrix
    m6  =  t(apply(hor_sym, 2, rev))
    m7  =  t(apply(m6, 2, rev))
    m8  =  t(apply(m7, 2, rev))
    
    # Combine all equivalent matrices into a list
    if(include_ == T){
      equivalences  =  list(land,as.vector(m2), as.vector(m3), as.vector(m4), as.vector(hor_sym), as.vector(m6), as.vector(m7), as.vector(m8))
    }else{
      equivalences  =  list(as.vector(m2), as.vector(m3), as.vector(m4), as.vector(hor_sym), as.vector(m6), as.vector(m7), as.vector(m8))
    }
    
    return(unique(equivalences))
  }
  
  generate_vectors_with_twos = function(N, x, y) {
    # Initialize an empty list to store the result vectors
    result_vectors <- list()
    
    # Special case: if there are no ones or no twos
    if (y == 0) {
      combinations_ones <- combn(N, x)
      for (i in 1:ncol(combinations_ones)) {
        vector <- rep(0, N)
        vector[combinations_ones[, i]] <- 1
        result_vectors[[length(result_vectors) + 1]] <- vector
      }
      return(result_vectors)
    }
    
    if (x - y == 0) {
      combinations_twos <- combn(N, y)
      for (i in 1:ncol(combinations_twos)) {
        vector <- rep(0, N)
        vector[combinations_twos[, i]] <- 2
        result_vectors[[length(result_vectors) + 1]] <- vector
      }
      return(result_vectors)
    }
    
    # General case: generate all combinations of positions for the 1s
    combinations_ones <- combn(N, x)
    
    # Iterate over each combination of positions for 1s
    for (i in 1:ncol(combinations_ones)) {
      # Create a base vector of zeros
      base_vector <- rep(0, N)
      
      # Place 1s in the specified positions
      base_vector[combinations_ones[, i]] <- 1
      
      # Generate all combinations of positions within the selected 1s to change to 2s
      positions_ones <- combinations_ones[, i]
      combinations_twos <- combn(positions_ones, y)
      
      # Iterate over each combination of positions for 2s
      for (j in 1:ncol(combinations_twos)) {
        # Create a copy of the base vector
        vector <- base_vector
        
        # Place 2s in the specified positions
        vector[combinations_twos[, j]] <- 2
        
        # Add the vector to the result list
        result_vectors[[length(result_vectors) + 1]] <- vector
      }
    }
    
    return(result_vectors)
  }
  
  landscape_to_hash = function(landscape) {
    return(digest(landscape, algo = 'xxhash32'))
  }
  
  low_lev_land4 = function(nb_nonz, nb2, size. = size) {
    #
    potential_landscapes = generate_vectors_with_twos(size.^2, nb_nonz, nb2)
    unique_landscapes = list()
    seen_hashes =  new.env(hash = TRUE, parent = emptyenv())
    
    for (landscape in potential_landscapes) {
      matrix_landscape = matrix(landscape, nrow = size.)
      equiv_landscapes = equi_landscape(matrix_landscape)
      
      is_unique = TRUE
      for (equiv in equiv_landscapes) {
        hash_equiv = landscape_to_hash(equiv)
        if (exists(hash_equiv, envir = seen_hashes)) {
          is_unique = FALSE
          break
        }
      }
      
      if (is_unique) {
        unique_landscapes[[length(unique_landscapes) + 1]] = landscape
        for (equiv in equiv_landscapes) {
          assign(landscape_to_hash(equiv), TRUE, envir = seen_hashes)
        }
      }
    }
    unique_landscapes = Matrix(do.call(rbind, unique_landscapes), sparse = T)
    writeMM(unique_landscapes,
            file = paste0(file_path_, 'data/landscape', size.,'x', size., '/landscape', size.,"_", nb_nonz,'_', nb2,'.mtx'))
    #return(unique_landscapes)
  }
  
  # Generate landscapes
  all_objects = ls(envir = .GlobalEnv)
  all_functions = all_objects[sapply(all_objects, function(x) is.function(get(x)))]
  
  
  candidates_full = list()
  
  for(nb_nonz_ in 1:size^2){
    nb2_ = seq(0,nb_nonz_)
    candidates_ = mapply(list, rep(nb_nonz_, nb_nonz_+1), nb2_, SIMPLIFY = F)
    candidates_full = append(candidates_full, candidates_)
  }
  
  start = Sys.time()
  
  cl = makeCluster(numCores)
  
  clusterExport(cl, varlist = c(all_functions, 'nb2_', 'nb_nonz_', "size", "file_path_"))
  clusterEvalQ(cl,{
    library(Matrix)
    library(digest)
  })
  
  result = parLapply(cl, candidates_full, function(inputs){
    low_lev_land4(inputs[[1]], inputs[[2]])
  })
  
  stopCluster(cl)
  
  
  #for(i in 1:length(result)){
  #  writeMM(result[i][[1]], file = paste0(file_path_, 'data/landscape', size,'x', size, '/landscape', size,"_", candidates_full[i][[1]][[1]],'_', candidates_full[i][[1]][[2]],'.mtx'))
  #}
  
  #for(nb2_ in 0:nb_nonz_){
  #  dat_ = low_lev_land4(nb_nonz_, nb2_)
  #  writeMM(dat_, file = here('landscapes_4x4', paste0("landscapes4_", nb_nonz_, '_', nb2_,'.mtx')))
  #}
  print(paste('Step', size, 'took', round(Sys.time() - start,3)))
}

