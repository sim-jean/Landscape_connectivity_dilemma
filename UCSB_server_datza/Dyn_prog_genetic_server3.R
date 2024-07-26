#### Dynamic programming for small scale landscapes ####
#### using fixed short term temporal horizons and   ####
#### solving heuristics (genetic algorithm)         ####
rm(list=ls())

# 0. Install packages ####

runner_ = 'UCSB'

required_packages = c('magrittr', 'dplyr', 'reshape2', 
                      'Matrix', 'digest', "parallel", "stringr","GA", 'stats')

if(runner_ == 'local'){
  library(here)
  working_path = paste0(here(),'/')
  numCores = 10
  
  required_packages = c('magrittr', 'dplyr', 'reshape2', 
                        'Matrix', 'digest', "parallel", "stringr","GA", 'stats')
  
  
  install_if_missing =function(package_){
    if (!package_ %in% installed.packages()) {
      install.packages(package_)
    }
    library(package_, character.only = TRUE)
  }
  
  sapply(required_packages, install_if_missing)
  
  
}else if(runner_ == 'INARI'){
  
  working_path = "/home/jean/connectivity_dilemma/data/"
  .libPaths("/home/jean/connectivity_dilemma/code/path_for_packages")
  options(repos = c(CRAN = "https://cran.r-project.org"))
  numCores = 10
  
  
  install_if_missing =function(package_){
    if (!package_ %in% installed.packages()) {
      install.packages(package_)
    }
    library(package_, character.only = TRUE)
  }
  
  sapply(required_packages, install_if_missing)
  
  
}else if(runner_ == 'UCSB'){
  working_path = '/home/simonjean/connectivity_2024/data/'
  numCores = as.numeric(Sys.getenv("SLURM_NTASKS"))
  user_lib = "~/R/library"
  
  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE)
  }
  
  .libPaths(user_lib)
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  
  install_if_missing <- function(package_) {
    if (!require(package_, character.only = TRUE)) {
      tryCatch({
        install.packages(package_, repos = "https://cloud.r-project.org/", lib = user_lib)
      }, error = function(e) {
        message(paste("Failed to install package:", package_, "\n", e))
      })
    }
    library(package_, character.only = TRUE, lib.loc = user_lib)
  }
  
  required_packages = c('RcppArmadillo', 'crayon', 'GA', 'foreach', 'iterators', required_packages)
  sapply(required_packages, install_if_missing)
  
  print('Set up on UCSB Pod finished!')
  
  sapply(required_packages, install_if_missing)
  
  print('Set up on UCSB Pod finished!')
}

# I. Parameters #### 
## a. Computational parameters
size = 3
max_budget = round(size^2/5)+1
budget = max_budget

#max_biodiv declared later
index1= size^2
index2 = size^2+1
index3 = 2*size^2
index4 = 2*size^2+1
index5 = 3*size^2
index6 = 3*size^2+1
index7 = 4*size^2
index8 = index7 +1
index9 = 5*size^2

popSize_ = 100
maxIter_ = 100

print('Primary parameters loaded')
print(paste('Population size is', popSize_))
print(paste('Number of iterations is', maxIter_))
## b. Meta parameters

set.seed(123)

data_path = paste0(working_path, 'landscape', size, 'x', size,'/landscapes/')
data_path2 = paste0(working_path, 'landscape', size, 'x', size,'/')


landscapes_to_treat = list.files(data_path)

# Order the files by size

#=print(landscapes_to_treat[1])
#landscapes_to_treat = landscapes_to_treat[1]
if(length(landscapes_to_treat)>0){
  print(paste0('Landscape files succesfully identified \n
        There are ', length(landscapes_to_treat), ' files to treat'))
  
}else{
  print("Error in identifying landscapes")
  print(paste(data_path, "may be inaccessible"))
}

sample_size = 2000
effective_sample = 0

print("Metaparameters loaded :")
print(paste0('Max sample size is: ', sample_size))
print(paste0('Number of cores is : ', numCores))



# II. Import functions #####
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

# III. Set up data ####
#biodivs_ = seq(0,1,.1)*max_biodiv



for(const_ in seq(.2,1,.1)*max_biodiv){
  for(file_name in landscapes_to_treat){
    
    landscape_ = readMM(paste0(data_path, file_name))
    
    print(paste(file_name, "data loaded"))
    
    #storage = paste0(data_path2, 'outputs/', strsplit(file_name, '.mtx')[[1]], '/')
    storage = paste0(data_path2, 'outputs/biod_', round(max_biodiv/const_),'/')
    
    if(!(file.exists(storage))){
      dir.create(storage)
    }
    
    storage2 = paste0(data_path2, 'outputs/utils/')
    if(!(file.exists(storage2))){
      dir.create(storage2)
    }
    
    print("Storage accessed or created")
    
    if(nrow(landscape_) > sample_size ){
      print('Too many landscapes : sampling')
      if(!file.exists(paste0(storage2, 'ids_',strsplit(file_name, '.mtx')[[1]],".csv"))){
        
        ids_ = generate_unique_uniform(sample_size, 1, nrow(landscape_))
        write.csv(ids_,paste0(storage2, 'ids_',strsplit(file_name, '.mtx')[[1]],".csv"), row.names= F)
        print('Initiated sampling and stored sample ID')
      }else{
        ids_ = read.csv(paste0(storage2, 'ids_',strsplit(file_name, '.mtx')[[1]],".csv"))
        print('Recovered sample ID')
      }
      ids_ = as.vector(unname(ids_))
      #print(ids_)
      landscape_ = landscape_[ids_, ,drop=F]
    }
    
    print(paste('For this run, we have', nrow(landscape_), 'runs'))
    effective_sample = effective_sample + nrow(landscape_)
    
    landscape_list = vector('list', nrow(landscape_))
    
    for(row_ in 1:nrow(landscape_)){
      landscape_list[[row_]] = landscape_[row_,]
    }
    rm(landscape_)
    print('Landscape data converted to list')
    # IV. Run parallel analysis. 
    
    # a. declare fitness function to be optimized by GA, eg taking treatments as inputs
    A7 = function(treatments_, 
                  landscape_,
                  B_ = budget,
                  constant_ = const_){
      # First try is landscape is 3*N^2
      # Generate indexes:
      
      
      # Generate treatments
      treatment0 = treatments_[1:index1]
      treatment1 = treatments_[index2:index3]
      treatment2 = treatments_[index4:index5]
      treatment3 = treatments_[index6:index7]
      treatment4 = treatments_[index8:index9]
      
      #Landscape and scores : 
      
      # Period 0 : 
      landscape_treated0 = landscape_*(1 - treatment0)
      score_fuel0 = mature_score(landscape_treated0, "F")
      score_biod0 = mature_score(landscape_treated0)
      
      # Period 1:
      suppressWarnings({
        landscape_treated1 = age_dyn(landscape_treated0)*(1 - treatment1)
      })
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
      # Period 4:
      suppressWarnings({
        landscape_treated4 = age_dyn(landscape_treated3)*(1 - treatment4)
      })
      score_fuel4 = mature_score(landscape_treated4, "F")
      score_biod4 = mature_score(landscape_treated4)
      
      # Constraints : 
      # Period0
      constraint_budget0 = B_ - sum(treatment0)
      constraint_biod0   = constant_ - score_biod0
      
      penalty1_0 = ifelse(constraint_budget0 >=0, 0, 10^4*abs(constraint_budget0))
      penalty2_0 = ifelse(constraint_biod0 <= 0, 0, 10^4*abs(constraint_biod0))
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
      # Period4
      constraint_budget4 = B_ - sum(treatment4)
      constraint_biod4   = constant_ - score_biod4
      
      penalty1_4 = ifelse(constraint_budget4 >=0, 0, 10^4*abs(constraint_budget4))
      penalty2_4 = ifelse(constraint_biod4 <= 0, 0, 10^4*abs(constraint_biod4))
      
      value_ = -(score_fuel0 + score_fuel1 + score_fuel2 + score_fuel3 + score_fuel4+
                   penalty1_0 + penalty1_1 + penalty1_2 + penalty1_3 +penalty1_4+
                   penalty2_0 + penalty2_1 + penalty2_2 + penalty2_3 +penalty2_4)
      
      return(as.numeric(value_))
    } 
    
    # b. declare function that takes landscape as input
    A7_ga = function(landscape_){
      
      A7_wrapper = function(treatments_){
        A7(treatments_, landscape_, B_ = budget, constant_= const_)
      }
      
      suggestions_ = matrix(rbinom(popSize_ * 5*size^2, 1, budget/size^2*0.9), nrow = popSize_, ncol = 5*size^2)
      
      #assign('landscape', landscape, envir = .GlobalEnv)
      
      ga_result = ga(
        type = "binary",
        fitness = A7_wrapper,
        nBits = 5*size^2,       # Number of bits in the binary vector (size of the landscape)
        maxiter = maxIter_,   # Number of generations
        popSize = popSize_,    # Population size
        pcrossover = 0.8,# Crossover probability
        pmutation = 0.2, # Mutation probability
        suggestions = suggestions_,
        monitor = F
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
    
    
    if(length(landscapes_to_treat)<numCores){
      numCores_ = length(landscapes_to_treat)
    }else{
      numCores_ = numCores
    }
    
    cl = makeCluster(numCores_)
    
    clusterExport(cl, varlist = c(all_objects[!(all_objects %in% c('landscapes_to_treat'))]))
    clusterEvalQ(cl,{
      library(Matrix)
      library(digest)
      library(magrittr)
      library(dplyr)
      library(GA)
    })
    
    
    
    
    # Parallel processing is not working and I dont get why. It worked perfectly a minute earlier
    result = parLapply(cl, landscape_list, A7_ga)
    
    stopCluster(cl)
    
    print(Sys.time() - start_)
    
    # Save data
    # first size^2 bits are the landscape, 3*size^2 last bits are the treatments. 
    result = t(Matrix(unlist(result), ncol = length(landscape_list), sparse=T))
    writeMM(result, paste0(storage,'results_', file_name ))
  }
  print(paste("Final analysis is done, effective sample size of ", effective_sample, "landscapes analyzed with biodiv ", const_))
}


for(const_ in seq(.1,1,1)*max_biodiv){
  for(file_name in landscapes_to_treat){
    
    landscape_ = readMM(paste0(data_path, file_name))
    
    print(paste(file_name, "data loaded"))
    
    #storage = paste0(data_path2, 'outputs/', strsplit(file_name, '.mtx')[[1]], '/')
    storage = paste0(data_path2, 'outputs/biod_', round(max_biodiv/const_),'/')
    
    if(!(file.exists(storage))){
      dir.create(storage)
    }
    
    storage2 = paste0(data_path2, 'outputs/utils/')
    if(!(file.exists(storage2))){
      dir.create(storage2)
    }
    
    print("Storage accessed or created")
    
    if(nrow(landscape_) > sample_size ){
      print('Too many landscapes : sampling')
      if(!file.exists(paste0(storage2, 'ids_',strsplit(file_name, '.mtx')[[1]],".csv"))){
        
        ids_ = generate_unique_uniform(sample_size, 1, nrow(landscape_))
        write.csv(ids_,paste0(storage2, 'ids_',strsplit(file_name, '.mtx')[[1]],".csv"), row.names= F)
        print('Initiated sampling and stored sample ID')
      }else{
        ids_ = read.csv(paste0(storage2, 'ids_',strsplit(file_name, '.mtx')[[1]],".csv"))
        print('Recovered sample ID')
      }
      ids_ = as.vector(unname(ids_))
      #print(ids_)
      landscape_ = landscape_[ids_, ,drop=F]
    }
    
    print(paste('For this run, we have', nrow(landscape_), 'runs'))
    effective_sample = effective_sample + nrow(landscape_)
    
    landscape_list = vector('list', nrow(landscape_))
    
    for(row_ in 1:nrow(landscape_)){
      landscape_list[[row_]] = landscape_[row_,]
    }
    rm(landscape_)
    print('Landscape data converted to list')
    # IV. Run parallel analysis. 
    
    # a. declare fitness function to be optimized by GA, eg taking treatments as inputs
    myopic_parallel = function(landscape, 
                               const_ = const_, 
                               budget_ = budget){
      
      landscapes_myopic_result = vector('list', 5)
      landscapes_myopic_result[[1]] = landscape
      treatments_ = vector('list', 5)
      
      risks = c()
      
      for(t in 2:6){
        landscapes_here = landscapes_myopic_result[[t-1]]
        const_ = biod  
        
        loc_optim = function(treatments_, 
                             landscape_,
                             B_ = budget,
                             constant_ = const_){
          # First try is landscape is 3*N^2
          # Generate indexes:
          
          # Generate treatme
          
          #Landscape and scores : 
          treatment0 = treatments_
          
          # Period 0 : 
          landscape_treated0 = landscape_*(1 - treatment0)
          score_fuel0 = mature_score(landscape_treated0, "F")
          score_biod0 = mature_score(landscape_treated0)
          
          
          # Constraints : 
          # Period0
          constraint_budget0 = B_ - sum(treatment0)
          constraint_biod0   = constant_ - score_biod0
          
          penalty1_0 = ifelse(constraint_budget0 >=0, 0, 10^4*abs(constraint_budget0))
          penalty2_0 = ifelse(constraint_biod0 <= 0, 0, 10^4*abs(constraint_biod0))
          
          
          value_ = -(score_fuel0+ 
                       penalty1_0 + 
                       penalty2_0)
          
          return(as.numeric(value_))
        } 
        
        # b. declare function that takes landscape as input
        loc_ga = function(landscape_){
          loc_wrapper = function(treatments_){
            loc_optim(treatments_, landscape_, B_ = budget, constant_= const_)
          }
          
          suggestions_ = matrix(rbinom(75* size^2, 1, budget/size^2*0.9), nrow = 75, ncol = size^2)
          
          #assign('landscape', landscape, envir = .GlobalEnv)
          
          ga_result = ga(
            type = "binary",
            fitness = loc_wrapper,
            nBits = size^2,       # Number of bits in the binary vector (size of the landscape)
            maxiter = 75,   # Number of generations
            popSize = 75,    # Population size
            pcrossover = 0.8,# Crossover probability
            pmutation = 0.2, # Mutation probability
            suggestions = suggestions_,
            monitor = F
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
          to_ret = c(treatment_)
          return(to_ret)
        }
        
        res_ =   loc_ga(landscapes_here)
        
        landscape_treated0 = landscapes_here*(1 - res_)
        risks = append(risks, mature_score(landscape_treated0, "F"))
        
        landscapes_myopic_result[[t]] = suppressWarnings({age_dyn(landscapes_here*(1-res_))})
        treatments_[[t-1]] = res_
      }
      return(c(landscape, unlist(treatments_)))
    }
    
    # C. set up objects to pass onto
    all_objects = ls(envir = .GlobalEnv)
    all_functions = all_objects[sapply(all_objects, function(x) is.function(get(x)))]
    
    start_ = Sys.time() 
    
    if(length(landscapes_to_treat)<numCores){
      numCores_ = length(landscapes_to_treat)
    }else{
      numCores_ = numCores
    }
    
    cl = makeCluster(numCores_)
    
    clusterExport(cl, varlist = c(all_functions, "landscape_list",'const_', 'biod', "popSize_", "maxIter_", "size", 
                                  "index1", "index2", "index3", "budget", "adjacency_" ))
    clusterEvalQ(cl,{
      library(Matrix)
      library(digest)
      library(magrittr)
      library(dplyr)
      library(GA)
    })
    
    treatments_myopic_loc = parLapply(cl, landscapes10, myopic_parallel)
    
    stopCluster(cl)
    
    treatments_myopic_loc_m = t(Matrix(unlist(treatments_myopic_loc), ncol = length(treatments_myopic_loc), nrow = 96))
    
    print(paste("Myopic with biodiveristy value", biod, "is done"))
    print(Sys.time() - start_biod_)
    writeMM(treatments_myopic_loc_m, here('landscape4x4', 'sample_analysis', paste0('result_myopic_biod_', biod, '.mtx')))
  }
  print(paste("Final analysis is done, effective sample size of ", effective_sample, "landscapes analyzed with biodiv ", const_))
}
