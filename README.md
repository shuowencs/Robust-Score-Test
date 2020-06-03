6/3/2020 Shuowen Chen and Hiroaki Kaido

This readme file describes the matlab code for project "Robust Score Test for Incomplete Models with Nuisance Parameters". 

The current main program is debugging.m. The program serves two purposes:
1. Check the performance of estimated nuisance parameter
2. Check the performance of the simulated test statstic (using either the true nuisance parameter or the first step estimand)

For readability, I categorize function files in accordance to the simulation procedures. Each file is commented to facilitate the understanding. 


# 1. Data Generation
We first generate data. Under the null, model is complete (checking the size property), selection is specified as 'Complete'. When model is not complete (checking the power property), use one of the three selection mechanisms (iid, non iid, LFP) to determine outcome among multiple equilibria. 

Function file that implements this: data_generation.m

# 2. Estimate the nuisance parameter
The test statistic requires a plugged-in estimand of nuisance parameter. The function file that implement this estimation (under the null that beta=0) is rmle.m. We consider three score-based algorithms: BHHH, BFGS and L-BFGS-M. For L-BFGS-M we use the code by Stephen Becker (https://github.com/stephenbeckr/L-BFGS-B-C). We prefer L-BFGS-M since it features better performance in our finite sample simulations. 

Auxiliary functions:
   counting.m: counts the number of occurrence of x (covariates) and s (outcomes)
   compute_z_all: computes analytical score vectors
Note: to implement L-BFGS-M, need to download and compile mexfiles from Stephen Becker's GitHub website. 
   
# 3. Constructing the sup test statistic
The function file is stat.m. In finite sample, we adjust the test statistic construction using Tikhonov regularization. 

# 4. Size Properties
Compute the actual size of the test statistic across different sample sizes. Data generation under the null

The main function that implements this is: actualsize.m

Auxiliary function:
crt.m: computes the critical value of the limiting distribution of the test statistic

# 5. Power Properties
Compute the power of the test statistic across different sample sizes. Data generation under the alternative

The main function that implements this is: powerfcn.m



