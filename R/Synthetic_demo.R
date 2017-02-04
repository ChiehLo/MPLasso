library(MASS)
library(Matrix)
library(matrixcalc)
library(glasso)
library(ROCR)
library(huge)
library(MCMCpack)
library(glmnet)
library(parallel)


source('~/Desktop/HMP_project/PLasso/R/syn_data_generate.R');
source('~/Desktop/HMP_project/PLasso/R/selection.R');
source('~/Desktop/HMP_project/PLasso/R/algorithm.R');



source('~/Desktop/HMP_project/SpiecEasi-master/R/normalization.R')


source('~/Desktop/HMP_project/my_R/evaluation.R');
source("~/Desktop/HMP_project/my_R/CalculateCoOccurrence.R")
source('~/Desktop/HMP_project/my_R/real_data.R')



run = 1;
eva = 1;
#load(file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/synthetic/cluster_100_200.RData");


if (run ==1){
  # Simulation Configuration
  OTUnum <- 50; # number of OTU
  simulations <- 2e2; # number of sample
  strength <- 0.2; # correlation strength
  numIt <- 2; # number of runs
  sampleLabel <- matrix(0, numIt, OTUnum*OTUnum);
  

  #PLasso Configuration
  prior <- 0.5;
  rmsePL = matrix(0, numIt);
  icovPL = matrix(0, numIt);
  rmsefPL = matrix(0, numIt);
  icovfPL = matrix(0, numIt);
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
    icovPL[i] <- evaluationPL$ICOV_l0;
    rmsefPL[i] <- evaluationPL$RMSE_F;
    icovfPL[i] <- evaluationPL$ICOV_F;
  }
}

if(eva==1){
  
  tempLabel <- split(t(sampleLabel), rep(1:nrow(sampleLabel), each = ncol(sampleLabel)))
  tempPred <- split(t(predPL), rep(1:nrow(predPL), each = ncol(predPL)))
  rocPL <- ROCnew(tempLabel, tempPred)
  perfPL <- rocPL$perf;
  aucPL <- rocPL$auc;
  aucstdPL <- rocPL$aucsd
  print("PLasso:")
  print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", mean(rmsePL), sd(rmsePL), mean(rmsefPL), sd(rmsefPL), aucPL, aucstdPL))
  #save(perfgp, file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/temp/test.RData")
  #load(file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/temp/test.RData")

  #save(sampleLabel, predcc, predgl, predgp, predm, predr, predre, preds, perfcc, perfgl, perfgp, perfmb, perfr, perfre, perfs, file = "/Users/chiehlo/Desktop/HMP_project/my_R/performance/synthetic/scale_free_100_200.RData")
}



