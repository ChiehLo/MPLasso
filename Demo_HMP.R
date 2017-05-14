#####################################################################
# Demo of the HMP datasets 
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
library(matrixStats)


#source files required
scriptDir <- dirname(sys.frame(1)$ofile)
source(paste(scriptDir, '/R/syn_data_generate.R', sep = ''));
source(paste(scriptDir, '/R/selection.R', sep = ''));
source(paste(scriptDir, '/R/algorithm.R', sep = ''));
source(paste(scriptDir, '/R/evaluation.R', sep = ''));
source(paste(scriptDir, '/R/calculate_co_occurrence.R', sep = ''))
source(paste(scriptDir, '/R/real_data.R', sep = ''))
source(paste(scriptDir, '/R/normalization.R', sep = ''))



# five body sites considered in this demo
bodySites <- c('Anterior_nares', 'Buccal_mucosa', 'Stool', 'Supragingival_plaque', 'Tongue_dorsum')
bodySites <- c('Anterior_nares')
# flag for evaluating different datasets
HMASM = TRUE;
HMMCP = TRUE;
HMQCP = TRUE;

# reproducibility evaluation
repro = TRUE;

if (HMASM == TRUE){
  for (i in 1:length(bodySites)){
    print("HMASM:")
    print(bodySites[i]);
    bodySite = bodySites[i];
    prePath = paste(scriptDir, "/data/HMASM/", sep = '');
    interactionFile = "/output/species_interaction.csv";
    interactionPath = paste(prePath, bodySite, interactionFile, sep= ''); 
    interactionMatrix <- interaction(interactionPath); # read interaction file 
    
    occuPath = paste(prePath, bodySite, "/output/association.csv", sep = '');
    NoAssociation <- occurance(occuPath); # read occurance file
    
    countPath = paste(prePath, bodySite, "/output/SpeciesCount.txt", sep = '');
    sample <- read_OTU(level = "species", countPath = countPath, 1-NoAssociation, interactionMatrix); # read OTU
    
    prior = 0;
    interactionFlag = TRUE;
    selectLambda <- selection(sample, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
    resultPL <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
    
    nInt <- 2;
    reproduceErrorArrayPL = matrix(0,nInt, 4); 
    if (repro == TRUE){
      for (k in 1:nInt){
        numSample <- nrow(sample$sample)
        subSample <- sample;
        sampleIndex <- sort(sample.int(numSample, round(0.5*numSample)))
        subSample$sample <- subSample$sample[sampleIndex,]
        lassoSub <- fraction_data_process(subSample$sample, 1 - NoAssociation, interactionMatrix)
        selectLambda <- selection(lassoSub, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
        lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
        resultPLSample <- algorithm_select(lassoSub, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
        reproduceErrorPL <- repro_eval(resultPL$adj, resultPLSample$adj)
        reproduceErrorArrayPL[k,] <- reproduceErrorPL;
       
      }
      meanPL <- colMeans(1-reproduceErrorArrayPL)
      stdPL <- colSds(1-reproduceErrorArrayPL)
      print("MPLasso HMASM Evaluation Results: 25% (std), 50% (std), 75% (std), 100% (std)")
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", meanPL[1], stdPL[1], meanPL[2], stdPL[2], meanPL[3], stdPL[3], meanPL[4], stdPL[4]))
    }
  }  
}

if (HMMCP == TRUE){
  for (i in 1:length(bodySites)){
    print("HMMCP:")
    print(bodySites[i]);
    bodySite = bodySites[i];
    prePath = paste(scriptDir, "/data/HMMCP/v35/", sep = '');
    occuPath = paste(prePath, bodySite, "/association.csv", sep = '');
    noAssociation <- occurance(occuPath);
  
    countPath = paste(prePath, bodySite, "/genus_count_filter.csv", sep = '');
    sample <- read_OTU(level = "genus", countPath = countPath, 1 - noAssociation);
    
    prior = 0;
    interactionFlag = FALSE;
    selectLambda <- selection(sample, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
    resultPLMM <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
    nInt <- 2;
    reproduceErrorArrayPL = matrix(0,nInt, 4);
    if (repro == TRUE){
      for (k in 1:nInt){
        numSample <- nrow(sample$sample)
        subSample <- sample;
        sampleIndex <- sort(sample.int(numSample, round(0.5*numSample)))
        subSample$sample <- subSample$sample[sampleIndex,]
        lassoSub <- count_data_process(subSample$sample, 1 - noAssociation)
        selectLambda <- selection(lassoSub, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
        lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
        resultPLSample <- algorithm_select(lassoSub, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
        reproduceErrorPL <- repro_eval(resultPLMM$adj, resultPLSample$adj)
        reproduceErrorArrayPL[k,] <- reproduceErrorPL;     
      }
      meanPL <- colMeans(1-reproduceErrorArrayPL)
      stdPL <- colSds(1-reproduceErrorArrayPL)
      print("MPLasso HMMCP Evaluation Results: 25% (std), 50% (std), 75% (std), 100% (std)")
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", meanPL[1], stdPL[1], meanPL[2], stdPL[2], meanPL[3], stdPL[3], meanPL[4], stdPL[4]))
    }
    
  }  
}

if (HMQCP == TRUE){
  for (i in 1:length(bodySites)){
    print("HMQCP:")
    print(bodySites[i]);
    bodySite = bodySites[i];
    prePath = paste(scriptDir, "/data/HMQCP/v35/", sep = '');
    occuPath = paste(prePath, bodySite, "/association.csv", sep = '');
    noAssociation <- occurance(occuPath);
    
    countPath = paste(prePath, bodySite, "/genus_count_filter.csv", sep = '');
    sample <- read_OTU(level = "genus", countPath = countPath, 1 - noAssociation);
    
    prior = 0;
    interactionFlag = FALSE;
    selectLambda <- selection(sample, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
    lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
    resultPLMQ <- algorithm_select(sample, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
    nInt <- 2;
    reproduceErrorArrayPL = matrix(0,nInt, 4);
    if (repro == TRUE){
      for (k in 1:nInt){
        numSample <- nrow(sample$sample)
        subSample <- sample;
        sampleIndex <- sort(sample.int(numSample, round(0.5*numSample)))
        subSample$sample <- subSample$sample[sampleIndex,]
        lassoSub <- count_data_process(subSample$sample, 1 - noAssociation)
        selectLambda <- selection(lassoSub, type="MPLasso", prior = prior, interactionFlag = interactionFlag);
        lambdaOptimal <- (which.min(selectLambda)-1)*0.01 + 0.01;
        resultPLSample <- algorithm_select(lassoSub, lambdaMin = lambdaOptimal, prior = prior, type = "MPLasso", interactionFlag = interactionFlag);
        reproduceErrorPL <- repro_eval(resultPLMQ$adj, resultPLSample$adj)
        reproduceErrorArrayPL[k,] <- reproduceErrorPL;

      }
      meanPL <- colMeans(1-reproduceErrorArrayPL)
      stdPL <- colSds(1-reproduceErrorArrayPL)
      print("MPLasso HMQCP Evaluation Results: 25% (std), 50% (std), 75% (std), 100% (std)")
      print(sprintf("%.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f)", meanPL[1], stdPL[1], meanPL[2], stdPL[2], meanPL[3], stdPL[3], meanPL[4], stdPL[4]))
    }
    
  }  
}
