### Attempt at genetic algorithms ####

install.packages("GA")
library(GA)
library(landscapeR)
library(raster)


# I. Parameters

size = 100
budget = round(size^2/5)

# II. Generate landscape ####
m  = matrix(0, size, size)
r  = raster(m, xmn = 0, xmx = size, ymn = 0, ymx = size)

num1 = round(runif(1, 0, 2*size))
size1 = round(runif(num1, 1, .8*size))
rr = makeClass(r, num1, size1, val= 1)

num2 = round(runif(1, 0, 2*size))
if(num2 == 0){
  size2 == 0
}else{
  size2 = round(runif(num2, 1, .8*size))
}
rr = makeClass(rr, num2, size2, val = 2)
plot(rr)

# Set landscape 
landscape = Matrix(rr@data@values, sparse = T)

landscape = Matrix(2, nrow = 100, ncol= 100)
# Set adjacency
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

adjacency_ = kings_graph(size)[1][[1]]

# Set biodiv
max_biodiv = as.numeric(Matrix(2, nrow = 1, ncol = size^2, sparse = T)%*%
  adjacency_ %*%
  Matrix(2, ncol = 1, nrow = size^2, sparse = T))

### Static optimization ####

fitness_function <- function(treatments,
                             B_ = budget) {
  # Inputs : 
    # Treatment is a binary vector of size N^2
    # Landscape is a vector as well 
  
  # Apply treatments to the landscape
  treated_landscape = landscape * (1 - treatment) 
  
  # Compute the first quadratic form
  score1 = as.numeric(t(Matrix(ifelse(treated_landscape == 2, 1, 0),
                    sparse = T)) %*% 
    adjacency_ %*% 
    Matrix(ifelse(treated_landscape == 2, 1, 0), sparse = T))
  
  score2 = as.numeric(t(Matrix(ifelse(treated_landscape > 0, 1, 0), 
                    sparse = T)) %*% 
    adjacency_ %*% 
    Matrix(ifelse(treated_landscape > 0, 1, 0), sparse = T))
  
  # Compute the objective function
  if (score2 < const_){
    if(sum(treatments) > B_ * 1.05){
      return(-(score1 + 10^8 +10^9))
    }else{
      return(-(score1 + 10^8)) 
    }
  } else {
    if(sum(treatments)>B_ *1.05){
      return(-(score1 + 10^9))
    }else{
    return(-score1)
    }
  }
}

repair_function = function(pop) {
  max_ones = budget
  for (i in 1:nrow(pop)) {
    num_ones = sum(pop[i,])
    if (num_ones > max_ones) {
      ones_indices = which(pop[i,] == 1)
      excess_ones = sample(ones_indices, num_ones - max_ones)
      pop[i, excess_ones] = 0
    }
  }
  return(pop)
}

ga_result <- ga(
  type = "binary",
  fitness = fitness_function,
  nBits = size^2,       # Number of bits in the binary vector (size of the landscape)
  maxiter = 100,   # Number of generations
  popSize = 200,    # Population size
  pcrossover = 0.2,# Crossover probability
  pmutation = 0.6, # Mutation probability
  elitism = 0.1,  # Elitism rate
  suggestions = matrix(rbinom(200 * 10000, 1, budget/size^2*0.95), nrow = 200, ncol = 10000)
)
### Illustration to see if it works
summary(ga_result)

treatment_ = ga_result@solution[15,]

a = rr
a = setValues(a, matrix(landscape, nrow = size))
plot(a)


b = rr
b = setValues(b, matrix(landscape*(1 - treatment_), nrow = size))
plot(b)


