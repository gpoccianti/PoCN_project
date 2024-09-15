library(igraph)
library(compiler)
library(deSolve)

#Wwe don't implement the noise
SWARM_r <- function(network, size, dt=5e-2, parms) {
  #parms : contains all the parameters except the ones linked to network topology (e.g. adjacency matrix)
  
  if(is.null(E(network)$weight)) stop("No edge weights specified.")
  
  # Get the connected components
  components <- components(network)
  component_sizes <- components$csize
  
  # Identify single-node components and exclude them
  non_single_nodes <- which(component_sizes[components$membership] > 1)
  if (length(non_single_nodes) == 0) stop("All nodes are single-node components.")
  
  network <- induced_subgraph(network, vids = non_single_nodes)  # Keep only non-single nodes
  parms <- parms[non_single_nodes,] #Keep only the corresponding parameters
  
  #prepare time series and parameters vectors
  S <- c()    #this is the network state, consisting of 3*N variables
  M <- 3      #number of variables for the model: 1-> phase, 2->x coordinate, 3->y coordinates
  idxs <- list()
  W <- get.adjacency(network, attr="weight", sparse=F)
  Nodes <- vcount(network)
  
  
  
  for(m in 1:M){
    idxs[[m]] <- (1:Nodes-1)*M + m
  }
  
  #initialization
  S[idxs[[1]]] <- runif(Nodes, 0, 2*pi)#initial phases
  S[idxs[[2]]] <- runif(Nodes, -1, 1) #inital x positions
  S[idxs[[3]]] <- runif(Nodes, -1, 1) #initial y positions
  
  k <-1 #iteration step (see numerical integration)
  
 
   
  edp <- ensure_distinct_points(S, idxs) #ensures no initial overlapping
  S <- edp$S
  
  #   MODEL
  #SWARMALATORS (2d):
  # also the variant with 1/N or 1/k_i as prefactor for the interaction
  model <- function(t, S, parms) {
    dS <- c()
    

   edp <- ensure_distinct_points(S, idxs) 
    S <- edp$S
    
    #calculate distances
  d_info <- calculate_distances(S,idxs)
  x_ij <- d_info$x_diff #x[i]-x[j] matrix
  y_ij <- d_info$y_diff #y[i]-y[j] matrix
  d <- d_info$distances #distance matrix
  
  #compute interactions  
    coupling <- parms$K/degree(network)
    
    D <- as.matrix(W/d) #adjacendy matrix weighted by inverse distances
    D_x <- as.matrix((parms$J*W*x_ij)/d) #adjacendy matrix weighted by inverse distances and d along x
    D_xx <- as.matrix((W*x_ij)/(d^2)) #adjacendy matrix weighted by squared inverse distances and d along x
    D_y <- as.matrix((parms$J*W*y_ij)/d) #adjacendy matrix weighted by inverse distances and  d along y
    D_yy <- as.matrix((W*y_ij)/(d^2)) #adjacendy matrix weighted by squared inverse distances and d along x
    
    dS[idxs[[1]]] <- parms$w  + coupling * ( cos(S[ idxs[[1]] ]) * (D %*% sin(S[ idxs[[1]] ])) - sin(S[ idxs[[1]] ]) * (D %*% cos(S[ idxs[[1]] ])) )

    # Spatial dynamics: x-coordinate
    dS[idxs[[2]]] <- parms$v_x - 1 / degree(network) * (D_x %*% as.matrix(parms$A/parms$J) + 
                                                          cos(S[idxs[[1]]]) * (D_x %*% cos(S[idxs[[1]]])) +
                                                          sin(S[idxs[[1]]]) * (D_x %*% sin(S[idxs[[1]]])) - 
                                                          D_xx %*% parms$B)
    
    # Spatial dynamics: y-coordinate
    dS[idxs[[3]]] <- parms$v_y - 1 / degree(network) * (D_y %*% as.matrix(parms$A/parms$J) + 
                                                          cos(S[idxs[[1]]]) * (D_y %*% cos(S[idxs[[1]]])) +
                                                          sin(S[idxs[[1]]]) * (D_y %*% sin(S[idxs[[1]]])) - 
                                                          D_yy %*% parms$B)

    
  #  if (any(abs(dS) > 1e6)) {
  #    warning("Extreme value detected in dS, potential numerical instability.")
  #  }
  
    cat(sprintf("\rIteration #%d at time %f", k, t))
    flush.console()
    k <<- k+1
    
    list(dS)
  }
  
  times <- seq(0, size*dt, by = dt)
  
  multi <- ode(y = S, times = times, func = model, parms = parms)
  #print(paste("DEBUG", dim(multi)))
  
  #transform in the MultiTS format
  #get only the variables of the state
  theta <- list()
  x <- list()
  y <- list()
  
  for(m in 1:Nodes){
    theta[[m]] <- multi[1:nrow(multi), idxs[[1]][m]+1]
    x[[m]] <- multi[1:nrow(multi), idxs[[2]][m]+1]
    y[[m]] <- multi[1:nrow(multi), idxs[[3]][m]+1]
  }
  
  return(list(th=theta,x=x,y=y))
}
SWARM <- cmpfun(SWARM_r)


####################
# OTHER FUNCTIONS
####################


normalize_vec <- function(v_centr){
  return((v_centr - min(v_centr))/( max(v_centr) - min(v_centr) ))
}

vec2pal <- function(v_centr, mypal){
  val_idxs <- 1 + floor((length(mypal)-1)*normalize_vec(v_centr))
  return(mypal[val_idxs])
}

plot.MultiTS <- function(MultiTS){
  library(ggplot2)
  
  dat <- data.frame()
  
  tseq <- 1:length(MultiTS[[1]])
  
  for(m in 1:length(MultiTS)){
    dat <- rbind(dat, data.frame(node=m, time=tseq, value=MultiTS[[m]]))
  }
  
  return( 
    ggplot(dat, aes(time, sin(value), group=node, color=node)) + theme_bw() + theme(panel.grid=element_blank()) +
      geom_line(alpha=0.3) + 
      scale_color_viridis_c() +
      xlab("Time") + ylab("sin(theta)")
  )
}

MultiTS2Matrix <- function(MultiTS){
  N <- length(MultiTS)
  M <- length(MultiTS[[1]])
  
  A <- matrix(0, nrow=N, ncol=M)
  
  for(n in 1:N){
    A[n,] <- MultiTS[[n]]
  }
  
  return(A)
}

#########


ensure_distinct_points <- function(S, idxs, min_distance = 1e-6, repulsion_strength = 0.01) {
  x_pos <- S[idxs[[2]]]
  y_pos <- S[idxs[[3]]]
  
  # Calculate pairwise distances
  x_diff <- outer(x_pos, x_pos, "-")
  y_diff <- outer(y_pos, y_pos, "-")
  dist_matrix <- sqrt(x_diff^2 + y_diff^2)
  
  # Find pairs that are too close
  too_close <- which(dist_matrix < min_distance & dist_matrix > 0, arr.ind = TRUE)
  
  while (nrow(too_close) > 0) {
    for (pair in 1:nrow(too_close)) {
      i <- too_close[pair, 1]
      j <- too_close[pair, 2]
      
      # Calculate the direction of repulsion
      repulsion_dir_x <- x_diff[i, j] / dist_matrix[i, j]
      repulsion_dir_y <- y_diff[i, j] / dist_matrix[i, j]
      
      # Apply repulsion to both nodes
      x_pos[i] <- x_pos[i] + repulsion_strength * repulsion_dir_x
      y_pos[i] <- y_pos[i] + repulsion_strength * repulsion_dir_y
      x_pos[j] <- x_pos[j] - repulsion_strength * repulsion_dir_x
      y_pos[j] <- y_pos[j] - repulsion_strength * repulsion_dir_y
    }
    
    # Recalculate distances
    x_diff <- outer(x_pos, x_pos, "-")
    y_diff <- outer(y_pos, y_pos, "-")
    dist_matrix <- sqrt(x_diff^2 + y_diff^2)
    
    # Update pairs that are too close
    too_close <- which(dist_matrix < min_distance & dist_matrix > 0, arr.ind = TRUE)
  }
  
  # Update S with new positions
  S[idxs[[2]]] <- x_pos
  S[idxs[[3]]] <- y_pos
  
  return(list(S = S, x_diff = x_diff, y_diff = y_diff))
}

calculate_distances <- function(S, idxs) {
  # Extract x and y positions based on indices
  x_pos <- S[idxs[[2]]]
  y_pos <- S[idxs[[3]]]
  
  # Calculate pairwise differences
  x_diff <- outer(x_pos, x_pos, "-")
  y_diff <- outer(y_pos, y_pos, "-")
  
  # Compute the pairwise distances
  dist_matrix <- sqrt(x_diff^2 + y_diff^2)
  
  # Ignore diagonal terms by setting them to Inf
  diag(dist_matrix) <- Inf
  
  # Return the difference matrices and distance matrix
  return(list(distances = dist_matrix, x_diff = x_diff, y_diff = y_diff))
}


