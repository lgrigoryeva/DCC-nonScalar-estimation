# DCC-nonScalar-estimation

This repository contains functions that carry out the maximum likelihood estimation of scalar and non-scalar DCC models. The nonlinear positive definiteness constraints are treated via with the Bregman divergences. All the computations are provided in the paper:

Bauwens, L., Grigoryeva, L. and Ortega, J.-P. [2015] Estimation and empirical performance of non-scalar dynamic conditional correlation models. To appear in  Computational Statistics and Data Analysis. doi:10.1016/j.csda.2015.02.013. http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2407652

- runDCCestimation.m allows to estimate scalar DCC, rank one deficient DCC, rank two deficient DCC, Hadamard DCC, Almon DCC, and Almon shuffle DCC models for the sample data file: exampleData.mat

REMARK: Notice that the parameters for the estimateDCC function have to be tuned. 
