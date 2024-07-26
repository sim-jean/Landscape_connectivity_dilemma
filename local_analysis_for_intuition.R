# Sample analysis locally #

library(stats)
library(dplyr)
library(magrittr)
library(Matrix)
library(GA)

working_path = paste0(here(),'/')
numCores = 10
size = 4
budget = 4
list_landscape = vector('list', 20)

for(i in 1:20){
  list_landscape[[i]] = runif(16, -0.45, 2.45)
}


kings_graph  =  function(n, timer_ = T) {
  if(timer_ == T){
    start = Sys.time()
  }
  # Initialize the adjacency matrix with zeros
  adj_matrix  =  Matrix(0, nrow = n*n, ncol = n*n, sparse = T)
  # Define potential moves (king moves in chess)
  moves  =  rbind(c(-1, -1), c(-1, 0), c(-1, 1), c(0, -1), c(0, 1), c(1, -1), c(1, 0), c(1, 1), c(0,0))
  
  # Helper function to safely add edges
  add_edge  =  function(i, j) {
    if (i >= 1 && i <= n*n && j >= 1 && j <= n*n) {
      adj_matrix[i,j]  =  1
      adj_matrix[j,i]  =  1
    }
    return(adj_matrix)
  }
  
  # Map each square to its neighbors
  for (row in 1:n) {
    for (col in 1:n) {
      index  =  (row - 1) * n + col
      
      for (move in 1:nrow(moves)) {
        neighbor_row  =  row + moves[move, 1]
        neighbor_col  =  col + moves[move, 2]
        if (neighbor_row > 0 && neighbor_row <= n && neighbor_col > 0 && neighbor_col <= n) {
          neighbor_index  =  (neighbor_row - 1) * n + neighbor_col
          if(adj_matrix[index,neighbor_index] == 0){
            adj_matrix = add_edge(index, neighbor_index) 
          }
        }
      }
    }
  }
  if(timer_ == T){
    end  =  Sys.time()
    time_taken  =  end - start
    return(list(matrix = adj_matrix, time = time_taken))
  } else {
    return(adj_matrix)
  }
}

if(file.exists(paste0(working_path,"adjacency_", size, ".mtx"))){
  adjacency_ = readMM(paste0(working_path,"adjacency_", size, ".mtx"))
}else{
  adjacency_ = kings_graph(size)[1][[1]]
  writeMM(adjacency_, file = paste0(working_path,"adjacency_", size, ".mtx"))
}

print('Adjacency loaded')

mature_cells = function(matrix){
  to_ret = Matrix(as.numeric(matrix>0), ncol = 1, sparse = T)
  return(to_ret)
}

burn_cells = function(matrix){
  to_ret = Matrix(as.numeric(matrix>1), ncol = 1, sparse = T)
  return(to_ret)
}

mature_score = function(matrix, 
                        type_ = 'B',
                        adj_ = adjacency_){
  if(type_ == 'B'){
    return(as.numeric(t(mature_cells(matrix)) %*% adj_ %*% mature_cells(matrix)))
  }else{
    return(as.numeric(t(burn_cells(matrix)) %*% adj_ %*% burn_cells(matrix)))
  }
}

age_dyn = function(landscape){
  landscape = Matrix(landscape, ncol = 1, sparse = T)
  #treatment = Matrix(treatment, ncol = 1, sparse = T)
  
  #augmented_mat = landscape + 1
  #capped_mat = pmin(augmented_mat, 2)
  
  #treatment_ = 1 - treatment
  #next_ = landscape * treatment_
  
  next_ = landscape + 1
  next_ = pmin(next_, 2)
  
  next_ = Matrix(as.numeric(next_), nrow = 1, byrow= T, sparse = T)
  return(next_)
}

# Functions for treatment allocation 

degrees = function(landscape){
  return(as.vector(t(Matrix(ifelse(landscape == 2, 1, 0),sparse = T)) %*% 
                     adjacency_))
  
}

generate_binary_matrix = function(N, length_of_vector, preferential_locations, X) {
  # N: number of binary vectors to generate
  # length_of_vector: the length of each binary vector
  # preferential_locations: indices where the 1s are preferentially located
  # X: maximum number of 1s in each vector
  
  # Initialize a matrix to store the generated binary vectors
  binary_matrix <- matrix(0, nrow = N, ncol = length_of_vector)
  
  for (i in 1:N) {
    # Initialize a binary vector with all 0s
    binary_vector <- rep(0, length_of_vector)
    
    # Randomly select the number of 1s to place, not exceeding X
    num_ones = X
    
    # Ensure preferential locations are not greater than the number of available positions
    num_preferential <- min(num_ones, length(preferential_locations))
    
    # Randomly select indices from preferential locations for the 1s
    if (num_preferential > 0) {
      selected_indices <- sample(preferential_locations, num_preferential)
    } else {
      selected_indices <- integer(0)
    }
    
    # Place 1s in the selected preferential locations
    binary_vector[selected_indices] <- 1
    
    # If needed, fill the remaining 1s in non-preferential locations
    remaining_ones <- num_ones - num_preferential
    if (remaining_ones > 0) {
      # Define non-preferential locations
      non_preferential_locations <- setdiff(1:length_of_vector, preferential_locations)
      
      # Randomly select remaining positions for the 1s
      remaining_indices <- sample(non_preferential_locations, remaining_ones)
      
      # Place the remaining 1s
      binary_vector[remaining_indices] <- 1
    }
    
    # Store the generated binary vector in the matrix
    binary_matrix[i, ] <- binary_vector
  }
  
  return(Matrix(binary_matrix, sparse = T))
}

max_biodiv = as.numeric(mature_score(Matrix(2,nrow = size^2)))

generate_unique_uniform <- function(n, min = 0, max = 1) {
  unique_draws <- numeric(0)
  while (length(unique_draws) < n) {
    new_draws <- round(runif(n - length(unique_draws), min, max))
    unique_draws <- unique(c(unique_draws, new_draws))
  }
  return(unique_draws[1:n])
}
print('Base functions imported')


# Analysis #

for(const_ in seq(0,1,.1)*max_biodiv){
  
  jvlivs = function(treatments_, 
                    landscape_,
                    B_ = budget,
                    constant_ = const_){
    # First try is landscape is 3*N^2
    # Generate indexes:
    
    
    # Generate treatments
    treatment1 = treatments_[1:index1]
    treatment2 = treatments_[index2:index3]
    treatment3 = treatments_[index4:index5]
    
    #Landscape and scores : 
    
    # Period 1 : 
    landscape_treated1 = landscape_*(1 - treatment1)
    score_fuel1 = mature_score(landscape_treated1, "F")
    score_biod1 = mature_score(landscape_treated1)
    
    # Period 2:
    suppressWarnings({
      landscape_treated2 = age_dyn(landscape_treated1)*(1 - treatment2)
    })
    score_fuel2 = mature_score(landscape_treated2, "F")
    score_biod2 = mature_score(landscape_treated2)
    
    # Period 3:
    suppressWarnings({
      landscape_treated3 = age_dyn(landscape_treated2)*(1 - treatment3)
    })
    score_fuel3 = mature_score(landscape_treated3, "F")
    score_biod3 = mature_score(landscape_treated3)
    
    
    # Constraints : 
    # Period1
    constraint_budget1 = B_ - sum(treatment1)
    constraint_biod1   = constant_ - score_biod1
    
    penalty1_1 = ifelse(constraint_budget1 >=0, 0, 10^4*abs(constraint_budget1))
    penalty2_1 = ifelse(constraint_biod1 <= 0, 0, 10^4*abs(constraint_biod1))
    # Period2
    constraint_budget2 = B_ - sum(treatment2)
    constraint_biod2   = constant_ - score_biod2
    
    penalty1_2 = ifelse(constraint_budget2 >=0, 0, 10^4*abs(constraint_budget2))
    penalty2_2 = ifelse(constraint_biod2 <= 0, 0, 10^4*abs(constraint_biod2))
    # Period3
    constraint_budget3 = B_ - sum(treatment3)
    constraint_biod3   = constant_ - score_biod3
    
    penalty1_3 = ifelse(constraint_budget3 >=0, 0, 10^4*abs(constraint_budget3))
    penalty2_3 = ifelse(constraint_biod3 <= 0, 0, 10^4*abs(constraint_biod3))
    
    value_ = -(score_fuel1 + score_fuel2 + score_fuel3 +
                 penalty1_1 + penalty1_2 + penalty1_3 +
                 penalty2_1 + penalty2_2 + penalty2_3)
    
    return(as.numeric(value_))
  } 
  # b. declare function that takes landscape as input
  jvlivs_ga = function(landscape_){
    
    jvlivs_wrapper = function(treatments_){
      jvlivs(treatments_, landscape_, B_ = budget, constant_= const_)
    }
    
    suggestions_ = matrix(rbinom(popSize_ * 3*size^2, 1, budget/size^2*0.9), nrow = popSize_, ncol = 3*size^2)
    
    #assign('landscape', landscape, envir = .GlobalEnv)
    
    ga_result = ga(
      type = "binary",
      fitness = jvlivs_wrapper,
      nBits = 3*size^2,       # Number of bits in the binary vector (size of the landscape)
      maxiter = maxIter_,   # Number of generations
      popSize = popSize_,    # Population size
      pcrossover = 0.8,# Crossover probability
      pmutation = 0.2, # Mutation probability
      suggestions = suggestions_,
      monitor = T
    )
    
    budg_eff = vector('list', nrow(ga_result@solution))
    
    for(i in 1:nrow(ga_result@solution)){
      budg_eff[i] = sum(ga_result@solution[i,])
    }
    
    budg_eff = unlist(budg_eff)
    treatment_id = which(budg_eff == min(budg_eff))
    
    if(length(treatment_id)>1){
      treatment_id = treatment_id[1]
    }
    treatment_ = unname(ga_result@solution[treatment_id,])
    #landscape
    to_ret = c(landscape_, treatment_)
    return(to_ret)
    
  }
  
  # C. set up objects to pass onto
  all_objects = ls(envir = .GlobalEnv)
  all_functions = all_objects[sapply(all_objects, function(x) is.function(get(x)))]
  
  start_ = Sys.time() 
  cl = makeCluster(numCores)
  
  clusterExport(cl, varlist = c(all_objects[!(all_objects %in% c('landscapes_to_treat'))]))
  clusterEvalQ(cl,{
    library(Matrix)
    library(digest)
    library(magrittr)
    library(dplyr)
    library(GA)
  })
  
  # Parallel processing is not working and I dont get why. It worked perfectly a minute earlier
  result = parLapply(cl, list_landscape, jvlivs_ga)
  
  stopCluster(cl)
  
  print(Sys.time() - start_)
}
