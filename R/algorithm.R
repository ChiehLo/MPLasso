#####################################################################
# PLasso implementation 
# You can implement your algorithm here
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

algorithm_select <- function(sample, lambdaMin = 0.01, lambdaMax = 10, prior = 1, type = "PLasso", priorType = TRUE, interactionFlag = FALSE){
  if(type=="PLasso"){
    result <- PLasso(sample, lambdaMin, lambdaMax, prior, priorType = priorType, interactionFlag = interactionFlag)
  }
  return(result);
}


PLasso <- function(sample, lambdaMin = 0.01, lambdaMax = 10, prior = 1, priorType, interactionFlag){
  s <- sample$var;
  adj <- sample$adj;
  if(interactionFlag==TRUE){
    interaction <- sample$inter;
  }
  p <- nrow(s);
  rho <- lambdaMin*matrix(1, p, p); #gLasso regularizer
  if (priorType == TRUE){
    rho <- lambdaMin*matrix(1, p, p); #gLasso regularizer
    for (i in 1:p){
      for (j in 1:p){
        if (abs(adj[i,j]) == 0 && runif(1, min=0, max=1) > prior){
          rho[i,j] <- lambdaMax;
          rho[j,i] <- lambdaMax;
        }
      }
    }
  }
  if(interactionFlag == TRUE){
    for (i in 1:p){
      for (j in 1:p){
        if (abs(interaction[i,j]) == 1 && runif(1, min=0, max=1) > prior && i!=j){
          rho[i,j] <- 0.00001;
          rho[j,i] <- 0.00001;
        }
      }
    }
  }
  a<-glasso(s, rho=rho)
  covOpt <- a$w;
  icovOpt <- a$wi;
  d <- 0;
  for(i in 1:p){
    d[i] <- 1/sqrt(covOpt[i,i]);
  }
  D <- diag(d);
  corrOpt <- D%*%covOpt%*%D;
  
  adjOpt <- matrix(0, p, p);
  for(i in 1:p){
    for (j in i:p){
      if (abs(icovOpt[i,j]) > 0){
        adjOpt[i,j] <- 1;
        adjOpt[j,i] <- 1;
      }
    }
  }
  
  
  return(list("cor" = corrOpt, "cov" = covOpt, "icov" = icovOpt, "adj" = adjOpt ));
}



