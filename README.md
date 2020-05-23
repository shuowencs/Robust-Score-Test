5/20/2020 Shuowen Chen and Hiroaki Kaido

This readme file describes the matlab code for project "Robust Score Test for Incomplete Models with Nuisance Parameters". 

The current main program is debugging.m. The program serves two purposes:
1. Check the performance of estimated nuisance parameter
2. Check the performance of the simulated test statstic (using either the true nuisance parameter or the first step estimand)

For readability, I categorize function files in accordance to the simulation procedures. Each file is commented to facilitate the understanding. 


# 1. Data Generation
We first generate data. Under the null, model is complete (checking the size property), selection is specified as 'Complete'. When model is not complete (checking the power property), use one of the three selection mechanisms (iid, non iid, LFP) to determine outcome among multiple equilibria. 

Function file that implements this: data_generation.m

# 2. Estimate the nuisance parameter
The test statistic requires a plugged-in estimand of nuisance parameter. The function file that implement this estimation (under the null that beta=0) is rmle.m. In small sample the estimand can be NA, and for these cases I use steepest ascent as an alternative. The code is in the rmle2.m. I commented their performances and issues in the debugging.m program. 

Auxiliary functions:
   counting.m: counts the number of occurrence of x (covariates) and s (outcomes)
   compute_z_all: computes analytical score vectors
   
# 3. Constructing test statistic
There are three functions that implements this. One is without regularization at all, while the other two feature different regularizations. 

(a) stat.m: without regularization
(b) tnDelta.m: apply regularization on the I_\delta\delta matrix
(c) tnVar.m: apply direct regularization on the orthogonalized variance covariance matrix

# 4. Size Properties
Compute the actual size of the test statistic across different sample sizes. Data generation under the null

The main function that implements this is: actualsize.m

Auxiliary function:
crt.m: computes the critical value of the limiting distribution of the test statistic

# 5. Power Properties
Compute the power of the test statistic across different sample sizes. Data generation under the alternative

The main function that implements this is: powerfcn.m



