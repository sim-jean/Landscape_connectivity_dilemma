########### Run sample analysis ###################
###################################################

required_packages = c('magrittr', 'dplyr', 'reshape2', 'here',
                      'Matrix', 'digest', "parallel", "stringr","GA", 'stats')

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

set.seed(123)

##### Parameters #####

size = 4

number_landscapes = 50
popSize_ = 200
maxIter_ = 250 

index1= size^2
index2 = size^2+1
index3 = 2*size^2
index4 = 2*size^2+1
index5 = 3*size^2
index6 = 3*size^2+1
index7 = 4*size^2
index8 = index7 +1
index9 = 5*size^2

max_budget = round(size^2/5)+1
budget = max_budget

list_landscape = vector('list', number_landscapes)

for(i in 1:number_landscapes){
  list_landscape[[i]] = round(runif(16, -0.45, 2.45))
}

#### Functions ######
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

degrees = function(landscape){
  return(as.vector(t(Matrix(ifelse(landscape == 2, 1, 0),sparse = T)) %*% 
                     adjacency_))
  
}

visualise_degrees = function(landscape, treatment){
  
  landscape_init = as.vector(landscape)
  landscape_post_treat = unname(landscape*(1 - treatment[1: index1]))
  landscape_init_degree = degrees(landscape_init)
  landscape_init_counter_degree = degrees(age_dyn(landscape_init))
  treatment_degree = as.vector(treatment[1:index1])*10
  
  landscape_init1 = as.vector(age_dyn(landscape_post_treat))
  landscape_post_treat1 = as.vector(landscape_init1*(1 - treatment[index2:index3]))
  landscape_init1_degree = degrees(landscape_init1)
  landscape_init1_counter_degree = degrees(age_dyn(landscape_init1))
  treatment1_degree = as.vector(treatment[index2:index3])*10
  
  landscape_init2 = degrees(as.vector(age_dyn(landscape_post_treat1)))
  landscape_post_treat2 = degrees(as.vector(landscape_init2*(1 - treatment[index4:index5])))
  landscape_init2_degree = degrees(landscape_init2)
  landscape_init2_counter_degree = degrees(age_dyn(landscape_init))
  treatment2_degree = as.vector(treatment[index4:index5])*10
  
  markers_1 = c(rep('initial', size^2),
                rep('counter0', size^2),
                rep('treat0', size^2),
                rep('initial1', size^2),
                rep('counter1', size^2),
                rep('treat1', size^2),
                rep('initial2', size^2),
                rep('counter2', size^2),
                rep('treat2', size^2))
  markers_2 = c(rep(0,3*size^2),
                rep(1,3*size^2),
                rep(2,3*size^2))
  
  dat_ =  as.data.frame(as.table(matrix(landscape_init_degree,size,size)))%>%
    mutate(Var1 = as.numeric(Var1),
           Var2 = as.numeric(Var2))
  
  for(x in c("landscape_init_counter_degree", 'treatment_degree', 'landscape_init1', "landscape_init1_counter_degree",
             'treatment1_degree', 'landscape_init2',"landscape_init2_counter_degree", "treatment2_degree")){
    
    candit = eval(parse(text = x))
    dat__ = as.data.frame(as.table(matrix(candit, size, size)))%>%
      mutate(Var1 = as.numeric(Var1),
             Var2 = as.numeric(Var2))
    dat_ = rbind(dat_, dat__)
  }
  
  dat_ = dat_%>%
    mutate(markers1 = factor(markers_1,
                             levels = c('initial','counter0', 'treat0', 'initial1','counter1',
                                        'treat1', 'initial2', 'counter2','treat2')),
           markers2 = factor(markers_2),
           Freq = as.factor(Freq))
  
  dat_ %>%
    ggplot(aes(x=Var1, y=Var2, fill = Freq))+
    geom_tile()+
    scale_fill_viridis_d()+
    facet_grid(markers2~markers1)+
    coord_fixed()+
    theme_minimal()
}

visualise_transition = function(landscape, treatment){
  # Evaluate 3 landscapes : 
  landscape_init = as.vector(landscape)
  landscape_post_treat = unname(landscape*(1 - treatment[1: index1]))
  landscape_init1 = as.vector(age_dyn(landscape_post_treat))
  landscape_post_treat1 = as.vector(landscape_init1*(1 - treatment[index2:index3]))
  landscape_init2 = as.vector(age_dyn(landscape_post_treat1))
  landscape_post_treat2 = as.vector(landscape_init2*(1 - treatment[index4:index5]))
  
  markers_1 = c(rep('initial', size^2),
                rep('post_treat0', size^2),
                rep('initial1', size^2),
                rep('post_treat1', size^2),
                rep('initial2', size^2),
                rep('post_treat2', size^2))
  markers_2 = c(rep(0,2*size^2),
                rep(1,2*size^2),
                rep(2,2*size^2))
  
  dat_ =  as.data.frame(as.table(matrix(landscape_init,size,size)))%>%
    mutate(Var1 = as.numeric(Var1),
           Var2 = as.numeric(Var2))
  
  for(x in c('landscape_post_treat', 'landscape_init1', 'landscape_post_treat1', 
             'landscape_init2', 'landscape_post_treat2')){
    
    candit = eval(parse(text = x))
    dat__ = as.data.frame(as.table(matrix(candit, size, size)))%>%
      mutate(Var1 = as.numeric(Var1),
             Var2 = as.numeric(Var2))
    dat_ = rbind(dat_, dat__)
  }
  
  dat_ = dat_%>%
    mutate(markers1 = factor(markers_1,
                             levels = c('initial', 'post_treat0', 'initial1',
                                        'post_treat1', 'initial2', 'post_treat2')),
           markers2 = factor(markers_2),
           Freq = as.factor(Freq))
  
  dat_ %>%
    ggplot(aes(x=Var1, y=Var2, fill = Freq))+
    geom_tile()+
    scale_fill_manual(values = green_palette)+
    facet_grid(markers2~markers1)+
    coord_fixed()+
    theme_minimal()
}

visualise_transition3 = function(landscape, treatment){
  # Evaluate 3 landscapes : 
  landscape_init = as.vector(landscape)
  landscape_post_treat = unname(landscape*(1 - treatment[1: index1]))
  landscape_init1 = as.vector(age_dyn(landscape_post_treat))
  landscape_post_treat1 = as.vector(landscape_init1*(1 - treatment[index2:index3]))
  landscape_init2 = as.vector(age_dyn(landscape_post_treat1))
  landscape_post_treat2 = as.vector(landscape_init2*(1 - treatment[index4:index5]))
  landscape_init3 = as.vector(age_dyn(landscape_post_treat2))
  landscape_post_treat3 = as.vector(landscape_init3*(1 - treatment[index6:index7]))
  landscape_init4 = as.vector(age_dyn(landscape_post_treat3))
  landscape_post_treat4 = as.vector(landscape_init4*(1 - treatment[index8:index9]))
  
  markers_1 = c(rep('initial', size^2),
                rep('post_treat0', size^2),
                rep('initial1', size^2),
                rep('post_treat1', size^2),
                rep('initial2', size^2),
                rep('post_treat2', size^2), 
                rep('initial3', size^2),
                rep('post_treat3', size^2),
                rep('initial4', size^2),
                rep('post_treat4', size^2))
  
  markers_2 = c(rep(0,2*size^2),
                rep(1,2*size^2),
                rep(2,2*size^2),
                rep(3,2*size^2),
                rep(4,2*size^2))
  
  dat_ =  as.data.frame(as.table(matrix(landscape_init,size,size)))%>%
    mutate(Var1 = as.numeric(Var1),
           Var2 = as.numeric(Var2))
  
  for(x in c('landscape_post_treat', 'landscape_init1', 'landscape_post_treat1', 
             'landscape_init2', 'landscape_post_treat2', "landscape_init3", 'landscape_post_treat3',
             "landscape_init4", 'landscape_post_treat4')){
    
    candit = eval(parse(text = x))
    dat__ = as.data.frame(as.table(matrix(candit, size, size)))%>%
      mutate(Var1 = as.numeric(Var1),
             Var2 = as.numeric(Var2))
    dat_ = rbind(dat_, dat__)
  }
  
  dat_ = dat_%>%
    mutate(markers1 = factor(markers_1,
                             levels = c('initial', 'post_treat0', 'initial1',
                                        'post_treat1', 'initial2', 'post_treat2', 'initial3', 'post_treat3', 'initial4', 'post_treat4')),
           markers2 = factor(markers_2),
           Freq = as.factor(Freq))
  
  dat_ %>%
    ggplot(aes(x=Var1, y=Var2, fill = Freq))+
    geom_tile()+
    scale_fill_manual(values = green_palette)+
    facet_grid(markers2~markers1)+
    coord_fixed()+
    theme_minimal()
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

#### Analysis ####

if(!(file.exists(here(working_path,"landscape4x4", 'sample_analysis', 'results_5peat_0.mtx')))){
  
  result_general = vector('list', length(seq(0,1,.1)))
  
  for(constraint_biod in seq(0,1,.1)*max_biodiv){
    
    const_ = constraint_biod
    
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
      landscape_treated0 = age_dyn(landscape_)*(1 - treatment0)
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
      
      suggestions_ = matrix(rbinom(popSize_ * 5*size^2, 1, budget/size^2*0.9), 
                            nrow = popSize_, ncol = 5*size^2)
      
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
    result = parLapply(cl, list_landscape, A7_ga)
    
    stopCluster(cl)
    
    print(Sys.time() - start_)
    result_general[[as.character(const_)]] = result
  }
  
  #    Save results
  
  for(const_ in as.character(seq(0,100,10))){
    loc_result = result_general[[const_]]
    to_save = t(Matrix(unlist(loc_result), ncol = number_landscapes, sparse = T))
    writeMM(to_save, here(working_path, 'landscape4x4', 
                          'sample_analysis', paste0('results_5peat_',const_, '.mtx')))
  }
  
}else{
  print("Results already computed")
}

if(!(file.exists(here(working_path,"landscape4x4", 'sample_analysis', 'result_myopic_biod_0.mtx')))){
  
  for(biod in seq(0,1,.1)*max_biodiv){
    
    start_biod_ = Sys.time()
    const_ = biod  
    
    myopic_parallel = function(landscape, 
                               const_ = const_, 
                               budget_ = budget){
      landscapes_myopic_result = vector('list', 8)
      landscapes_myopic_result[[1]] = landscape
      treatments_ = vector('list', 8)
      
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
          landscape_treated0 = age_dyn(landscape_)*(1 - treatment0)
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
    
    all_objects = ls(envir = .GlobalEnv)
    all_functions = all_objects[sapply(all_objects, function(x) is.function(get(x)))]
    
    numCores = 10
    
    cl = makeCluster(numCores)
    
    clusterExport(cl, varlist = c(all_functions, 'landscapes10','const_', 'biod', "popSize_", "maxIter_", "size", 
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
    
    treatments_myopic_loc_m = t(Matrix(unlist(treatments_myopic_loc), ncol = length(treatments_myopic_loc), nrow = 7*(size^2)))
    
    print(paste("Biodiveristy value", biod, "is done"))
    print(Sys.time() - start_biod_)
    writeMM(treatments_myopic_loc_m, here('landscape4x4', 'sample_analysis', paste0('result_myopic_biod_', biod, '.mtx')))
  }
}
