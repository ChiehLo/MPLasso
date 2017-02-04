#####################################################################
# Demo of the synthetic experiment
#
# @author Chieh Lo
# @date 01/26/2017
#####################################################################

# library required
library(MASS)
library(Matrix)
library(matrixcalc)
library(glasso)
library(ROCR)
library(huge)
library(MCMCpack)
library(glmnet)
library(parallel)

#source files required
scriptDir <- dirname(sys.frame(1)$ofile)
source(paste(scriptDir, '/R/syn_data_generate.R', sep = ''));
source(paste(scriptDir, '/R/selection.R', sep = ''));
source(paste(scriptDir, '/R/algorithm.R', sep = ''));
source(paste(scriptDir, '/R/evaluation.R', sep = ''));
source(paste(scriptDir, '/R/calculate_co_occurrence.R', sep = ''))
source(paste(scriptDir, '/R/normalization.R', sep = ''))


# flag for running and evaluation
run = TRUE;
eva = TRUE;



if (run ==TRUE){
  # Simulation Configuration
  OTUnum <- 50; # number of OTU
  simulations <- 2e2; # number of sample
  strength <- 0.2; # correlation strength
  numIt <- 2; # number of runs
  sampleLabel <- matrix(0, numIt, OTUnum*OTUnum);
  

  #PLasso Configuration
  prior <- 0.5;
  rmsePL = matrix(0, numIt);
  predPL = matrix(0, numIt, OTUnum*OTUnum);

  for (i in 1:numIt){
    
    #sample generation
    sample <- graph_select(OTUnum = OTUnum, Simulations = simulations, Strength = strength, type = "random");
    sampleLabel[i,] <- abs(as.vector(sample$adj));
    
    # model selection
    selectLambda <- selection(sample, type="PLasso", prior = prior);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;

    # optimal model
    resultPL <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "PLasso");
    predPL[i,] <- abs(as.vector(resultPL$icov));
    
    # calculate L1 distance
    evaluationPL <- evaluation(sample, resultPL);
    rmsePL[i] <- evaluationPL$RMSE_l0;
  }
}

if(eva==TRUE){
  tempLabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  tempPred <- split(t(predPL), rep(1:nrow(predPL), each = ncol(predPL)))
  rocPL <- ROC_eval(tempLabel, tempPred)
  perfPL <- rocPL$perf;
  aucPL <- rocPL$auc;
  aucstdPL <- rocPL$aucsd
  acc <- acc_eval(tempLabel, tempPred);
  print("PLasso Evaluation Results: rmse (std), acc (std), auc (std)")
  print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(rmsePL), sd(rmsePL), mean(acc), sd(acc), aucPL, aucstdPL))
}



