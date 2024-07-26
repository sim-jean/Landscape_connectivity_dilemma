#### Dynamic programming for server use ####

# 0. Install packages ####
required_packages = c('magrittr', 'dplyr', 'reshape2', 
                      'Matrix', 'digest', "parallel", "stringr")
  
install_if_missing =function(package_){
  if (!package_ %in% installed.packages()) {
    install.packages(package_)
  }
  library(package_, character.only = TRUE)
}

sapply(required_packages, install_if_missing)

# I. Declare key parameters, integral and meta parameters ####

size = 3

time_horizon = 10
const_ = 40

numCores = 10

budget = round(size^2/5)+1

# II. Declare functions ####

kings_graph  =  function(n, timer_ = T) {
  if(timer_ == T){
    start = Sys.time()
  }
  # Initialize the adjacency matrix with zeros
  adj_matrix  =  Matrix(0, nrow = n*n, ncol = n*n)
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

adjacency_ = kings_graph(size)$matrix
max_biodiv = as.numeric(Matrix(1, nrow = 1, ncol = size^2, sparse = T)%*%
                          adjacency_ %*%
                          Matrix(1, ncol = 1, nrow = size^2, sparse = T))


mature_cells = function(matrix){
  to_ret = Matrix(as.numeric(matrix<2), ncol = 1, sparse = T)
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
    return(t(mature_cells(matrix)) %*% adj_ %*% mature_cells(matrix))
  }else{
    return(t(burn_cells(matrix)) %*% adj_ %*% burn_cells(matrix))
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

to_optimize = function(landscape, constraint_){
  
  y = mature_score(landscape, type_='F') + ifelse(mature_score(landscape) <= constraint_, 10^5, 0)
  
  return(y)
}

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
    equivalences  =  list(as.vector(land), as.vector(m2), as.vector(m3), as.vector(m4), as.vector(hor_sym), as.vector(m6), as.vector(m7), as.vector(m8))
  }else{
    equivalences  =  list(as.vector(m2), as.vector(m3), as.vector(m4), as.vector(hor_sym), as.vector(m6), as.vector(m7), as.vector(m8))
  }
  
  return(unique(equivalences))
}

# II. Generate treatments ####

number_treats   = paste0('choose(', size^2 ,',', seq(1:budget),')' )
number_treats   = eval(parse(text= paste(number_treats, collapse = '+')))
potential_treat = vector('list', number_treats)


generate_combinations  =  function(N, k_max) {
  all_combinations  =  list()
  total_columns  =  0
  
  for (k in 1:k_max) {
    comb  =  combn(N, k)
    total_columns  =  total_columns + ncol(comb)
    binary_vectors  =  apply(comb, 2, function(indices) {
      vec  =  rep(0, N)
      vec[indices]  =  1
      return(vec)
    })
    all_combinations[[k]]  =  binary_vectors
  }
  
  # Combine all combinations into a single matrix
  combined_matrix  =  do.call(cbind, all_combinations)
  
  # Convert to sparse matrix
  sparse_comb_matrix  =  Matrix(combined_matrix, sparse = TRUE)
  
  return(sparse_comb_matrix)
}

treat = generate_combinations(size^2, budget)

# III. Set up storage space ####

suppressWarnings({
  # 1. Set time horizon
  
  # 2. Discretize search space : well, it is already discrete, it is the real of possible landscapes ####
  
  #file_path = here(paste0('landscapes_', size, 'x', size),'data', 'landscapes')
  file_path = paste0('/home/jean/connectivity_dilemma/data/landscape',size, 'x', size, '/landscapes_/')
  
  files_ = list.files(file_path)
  
  landscapes_ = vector('list', length(files_))
  
  for(i in 1:length(files_)){
    landscapes_[i] = readMM(paste0(file_path, '/', files_[i]))
  }
  
  extract_pairs <- function(x) {
    matches <- str_extract_all(x, "\\d+")[[1]]
    if (length(matches) >= 2) {
      return(c(matches[2:3]))
    } else {
      return(NA)
    }
  }
  
  files_id_ = lapply(files_, extract_pairs)
  files_id_final = vector('list', length= length(files_id_))
  for(i in 1:length(files_id_)){
    files_id_final[[i]] = paste0(as.character(unlist(files_id_[i])), collapse = '')
  }
  names(landscapes_) = files_id_final
  
  # 3. Value function store for the N possible stock values ####
  # Values are stored in lists, representing the number of non zero and twos
  # In each list, there are as many values as unique landscapes
  V_end = vector('list', length(landscapes_))
  names(V_end) = files_id_final
  
  for(i in 1:length(V_end)){
    loc_nrow = nrow(landscapes_[[i]])
    V_end[i] = Matrix(0, nrow = loc_nrow)
  }
  
  treats_ = vector('list', length(landscapes_))
  names(treats_) = files_id_final
  
  for(i in 1:length(treats_)){
    loc_nrow = nrow(landscapes_[[i]])
    treats_[i] = Matrix(0, nrow = loc_nrow)
  }
  
  # 4. Loop backwards over time, going from T to 1
  
  payoff = function(landscape, 
                    treatment, 
                    constraint_ = constraint){
    # Inputs : landscape, in vector form
    #          treatment, in vector form
    #          constraint_, numeric
    
    candidate_ = landscape*(1 - treatment)
    
    # Compute current value
    current_ = to_optimize(candidate_, constraint_)
    
    # Compute next period landscape
    x_next = age_dyn(candidate_)
    
    # Find where to pick the continuation value :
    # a. Find what matrix to look into
    nb2_next = sum(x_next==2)
    nonz_next = sum(x_next>0)
    lookup = paste0(as.character(nonz_next),as.character(nb2_next))
    
    # b. Once we have found the right matrix to look into, we can find which matrix is equivalent to the considered matrix
    # bb. May need to figure a way to have smarter lookout : take the big matrix and cut it. We'll see for 4x4
    to_check = landscapes_[[lookup]]
    
    # Generate equiv_landscapes
    equiv_land = equi_landscape(matrix(x_next, nrow=size), include_ = T)
    
    if(nrow(to_check) == 1){
      index_to_take=1
    }
    # Set the hashes
    seen_hashes =  new.env(hash = TRUE, parent = emptyenv())
    
    landscape_to_hash = function(land) {
      return(digest(land, algo = 'xxhash32'))
    }
    
    add_hash = function(land){
      assign(landscape_to_hash(land), TRUE, envir = seen_hashes)
    }
    # Compute all the hashes of 
    
    apply(to_check, 1, add_hash)
    
    # for each potential equiv, check if the hash is in the environment of hashes
    # When it is found, we have the right index to take the continuation value from
    for(equiv in equiv_land){
      hash_equiv = landscape_to_hash(equiv)
      if(exists(hash_equiv, envir = seen_hashes)) {
        index_to_take = which(ls(seen_hashes)==hash_equiv)
        break
      }
    }
    
    continuation = V_end[[lookup]][index_to_take]
    
    return(as.numeric(current_ + continuation))
  }
  
  # Now, we will need to look for every element of the landscapes, and find the optimal value for each of the
})

# IV. Initiate storage for different time horizons ####

Vs = vector('list', time_horizon+1)
Vs[[time_horizon+1]] = V_end

treatments_ = vector('list', time_horizon+1)

# V. Dynamic programming

dyn_prog_indiv = function(loc_landscapes){
  # Required : 
  # Landscapes_[i] : input space
  # treat : to be passed on
  # V : to be passed on
  
  
  namer_ = names(loc_landscapes)
  
  #loc_landscapes = loc_landscapes[[1]]
  # Initiate output
  V_out = Matrix(1, nrow = nrow(loc_landscapes), sparse = T)
  treat_out = Matrix(1, nrow = nrow(loc_landscapes), sparse = T)
  
  for(j in 1:nrow(loc_landscapes)){
    #nrow(loc_landscapes)){
    # Loop over each landscape
    
    loc_land = as.numeric(loc_landscapes[j,])
    payoffs_ = vector('list', ncol(treat))
    
    for(tr in 1:ncol(treat)){
      payoffs_[tr] = payoff(loc_land, treat[,tr], constraint_ = const_)
    }
    
    # Find where the minimum payoff is
    payoffs_ = unlist(payoffs_)
    value_attributed = unique(min(payoffs_))
    potential_treatment = which(payoffs_==value_attributed)
    
    dat_here = as.data.frame(matrix(treat[,potential_treatment], nrow=size^2))%>%
      colSums() # pick minimum budget
    names(dat_here) = potential_treatment
    chosen_treatment  = as.numeric(names(which(dat_here==min(dat_here))))
    
    if(length(chosen_treatment)>1){
      chosen_treatment = chosen_treatment[round(runif(1,1,length(chosen_treatment)))]
    }
    # Store the chosen treatment
    # Update Vend
    V_out[j,] = value_attributed
    treat_out[j,] = chosen_treatment
  }
  return(list(V_out, treat_out, namer_))
}

all_objects = ls(envir = .GlobalEnv)
all_functions = all_objects[sapply(all_objects, function(x) is.function(get(x)))]


for(t in time_horizon:1){
  print(paste('Period ', t, 'is starting'))
  
  # Initiate clusters
  start_ = Sys.time() 
  cl = makeCluster(numCores)
  
  clusterExport(cl, varlist = c(all_functions, 'landscapes_', 'V_end', 
                                "size", "budget", "const_", "treat", "adjacency_"))
  clusterEvalQ(cl,{
    library(Matrix)
    library(digest)
    library(magrittr)
    library(dplyr)
  })
  
  
  result = parLapply(cl, landscapes_, dyn_prog_indiv)
  stopCluster(cl)
  
  print(Sys.time() - start_)
  
  # Assign values
  
  for(name_ in names(result)){
    V_end[[name_]] = result[[name_]][[1]]
    treats_[[name_]]= result[[name]][[2]]
  }
  
  Vs[[t]] = V_end
  treatments_[[t]] = treats_
}

saveRDS(Vs, file = paste0('/home/jean/connectivity_dilemma/data/landscape',size, 'x', size, '/values_/values_biodiv', const_, '.rds'))
saveRDS(Vs, file = paste0('/home/jean/connectivity_dilemma/data/landscape',size, 'x', size, '/treatments_/treatments_biodiv', const_, '.rds'))



