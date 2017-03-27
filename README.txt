# MPLasso
################################################################################
File: README.md
A brief introduction about the usage for MPLasso
-------------------------------------------------------------------------------
Author: Chieh Lo (Carnegie Mellon University)
Email : chiehl@andrew.cmu.edu
Date  : 2017-03-22
-------------------------------------------------------------------------------
Files: 
               ./Demo_Synthetic.R: A demo for synthetic datasets. 
               ./Demo_HMP.R: A demo for HMP datasets.
               ./R/algorithm.R: Implementation of the MPLasso
               ./R/calculate_co_occurrence.R: Calculate the prior information based on microbial co-occurrence.
               ./R/evaluation.R: Evaluate the performacne of the MPLasso.
               ./R/normalization.R: Log-ratio transformation for microbial data. 
               ./R/real_data.R: HMP data preprocessing. 
               ./R/selection.R: Model selection for MPLasso based on BIC.
               ./R/syn_data_generate.R: Generation of synthetic data of five different graphs using huge package.
HMP datasets:
               ./data/HMASM/: Post-processing count data and the prior information obtained using @MInter. 
               ./data/HMMCP/: Post-processing count data and the prior information obtained using @MInter
               ./data/HMQCP/: Post-processing count data and the prior information obtained using @MInter
               For HMASM, it includes both species interactions and co-occurrance prior information.
-------------------------------------------------------------------------------

Quick start:
source Demo_Synthetic.R or Demo_HMP.R in Rstudio.

