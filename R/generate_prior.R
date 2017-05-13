generatePrior <- function(sample, prior, precision){
  p <- ncol(sample$sample);
  n <- nrow(sample$sample);
  adj <- sample$adj;
  rho_inter <- matrix(1, p, p);
  diag(rho_inter) <- 0
  rho_occur <- matrix(0, p, p);
  # type 1: know the interaction, i.e., know adj[i,j] == 1 
  if (prior >= 0 && precision > 0){
    precision_ratio <- precision
    # 1 for 
    rho_inter <- matrix(1, p, p); #gLasso regularizer
    diag(rho_inter) <- 0
    # generate the constriant matrix 
    count_penalized <- 0
    total_count <- 0
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 1 && runif(1, min=0, max=1) > prior && i != j){
          rho_inter[i,j] <- 0;
          #rho_inter[j,i] <- 0;
          if(runif(1, min=0, max=1) > (1- precision_ratio) ){
            count_penalized <- count_penalized + 1
            rho_inter[i,j] <- 1;
            #rho_inter[j,i] <- 1;
          }
          total_count <- total_count + 1
        }
      }
    }
    print(count_penalized)
    print(total_count)
    while (count_penalized > 0){
      row = round(runif(1, min = 1, max = p))
      col = round(runif(1, min = 1, max = p))
      if(row!=col && rho_inter[row, col]==1){
        rho_inter[row,col] = 0;
        #rho_inter[col,row] = 0;
        count_penalized <- count_penalized -1
      }
    }
    print(count_penalized)
  }
  
  # type 2: know non interacting edges, i.e., adj[i,j] == 0
  if (prior >= 0 && precision > 0){
    precision_ratio <- precision
    rho_occur <- matrix(0, p, p); #gLasso regularizer
    # generate the constriant matrix
    count_penalized <- 0
    total_count <- 0
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior && i != j){
          rho_occur[i,j] <- 1;
          #rho_occur[j,i] <- 1;
          if(runif(1, min=0, max=1) > (1- precision_ratio) ){
            count_penalized <- count_penalized + 1
            rho_occur[i,j] <- 0;
            #rho_occur[j,i] <- 0;
          }
          total_count <- total_count + 1
        }
      }
    }
    print('type2:')
    print(total_count)
    print(count_penalized)
    counts <- 5000000;
    while(count_penalized > 1000){
      row = round(runif(1, min = 1, max = p))
      col = round(runif(1, min = 1, max = p))
      if(row!=col && rho_occur[row, col]==0){
        rho_occur[row,col] = 1;
        #rho_occur[col,row] = 1;
        count_penalized <- count_penalized -1
      }
      counts <- counts - 1;
    }
    print(count_penalized)
  }
  
  if(prior>=0 && precision == 0){
    rho_inter <- matrix(1, p, p); #gLasso regularizer
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 1 && runif(1, min=0, max=1) > prior && i!=j){
          rho_inter[i,j] <- 0;
          #rho_inter[j,i] <- 0;
        }
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior && i != j){
          rho_occur[i,j] <- 1;
          #rho_occur[j,i] <- 1;
        }
      }
    }
  }
  return(list("rho_inter" = rho_inter, "rho_occur" = rho_occur));
}