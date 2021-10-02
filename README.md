10/1/2021 Shuowen Chen and Hiroaki Kaido

This readme file describes the matlab code for project "Robust Tests for Model Incompleteness with the Presence of Nuisance Parameters". 

# Monte Carlo

The main program for Monte Carlo simulation is monte_carlo.m. For readability, we categorize function files in accordance to the simulation procedures. Each file is commented to facilitate the understanding. 

## 1. Data Generation
We first generate data. Under the null, model is complete (checking the size property), selection is specified as 'Complete'. Under the alternative, the model is incomplete (checking the power property), use one of the three selection mechanisms (iid, non iid, LFP) to determine outcome data among multiple equilibria. 

Function file that implements this: data_generation.m

## 2. Estimate the nuisance parameter
The test statistic requires a plugged-in estimator of nuisance parameter. We consider three score-based algorithms: BHHH, BFGS and L-BFGS-M. For L-BFGS-M we use the code by Stephen Becker (https://github.com/stephenbeckr/L-BFGS-B-C). We prefer L-BFGS-M since it features better performance in our finite sample simulations. However, when sample size is large, the three algorithms provide similar estimation results. 

The function file that implements this estimation (under the null that beta=0): rmle.m. 

Auxiliary functions:
   counting.m: counts the number of occurrence of x (covariates) and s (outcomes)
   compute_z_all: computes analytical score vectors
   
Note: to implement L-BFGS-M for the first step estimation, need to download from Stephen Becker's GitHub website compile mex files in Matlab.
   
## 3. Constructing the test statistic
The function file is stat_ab.m. For finite-sample property, we regularize some components of the test statistic using suggestions from Andrews and Barwick (2012). 

## 4. Critical values

The function file that simulates critical values: crt_ab.m

$$ 5. Post model selection inference on nuisance parameters
Compare the test statistic to shrinkage critical value and decide if Wald test or robust test (e.g., Kaido, Molinari and Stoye, 2019) should be implemented. 

Functions that implement respective inference: get_Wald.m and get_KMS.m in the KMS folder. 

# Comparison with Bugni, Canay and Shi (2017)
We compare the robust score test with MR test by BCS (QE, 2017). We modify their online code to incorporate covariates. 

The main file is Simulation_BCS.m. Auxiliary functions are

[Modified]
(1) dataOp.m: computes standardized sample moment equalities and inequalities 
(2) Qn_function_X.m: defines sample criterion function according to MMM in Eq. (4.2) of BCS with covariates
(3) Qn_MR_function_X.m: sample criterion function for DR and PR inference 

[Hardcopy]
(4) S_function.m: called by Qn_function_X.m and Qn_MR_function_X to compute the S function in the main text of BCS
(5) phi_function.m: Defines the GMS penalization function for DR inference as in Eq. (2.11) of BCS

# Application

The main program for application is application.m. The data we use is from Kline and Tamer (2016, QE), which is available from QE website. 

We consider two alternatives: with and without discretization

In addition to crt_ab.m, other auxiliary functions are 
(1) rmleapp.m: estimates nuisance parameters under the null; 
(2) lhsscale.m: implements Latin hypercube sampling to generate multiple starting points for the estimation of nuisance parameters. 
(3) statapp.m: computes test statistic 
(4) rmle_nodis.m: rmle without discretization
(5) test_stat_no_discretization: compute test statistic without discretization

