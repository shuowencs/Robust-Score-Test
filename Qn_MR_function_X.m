%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   5/1/2021 Shuowen Chen and Hiroaki Kaido
%   This is taken from code by Bugni, Canay and Shi (QE) with some changes
%   (1) We have two type of moment conditions
%   (2) p, the number of moment inequalities, are now in the fcn body
%   (3) k, total number of moments, are now in the fcn body
%   (3) we use dataOp.m to produce studendized moments
%   (4) instead of naming theta_to_min and theta_H0, we directly input delta
%       and beta
%   Defines the sample criterion function for DR or PR inference
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = Qn_MR_function_X(delta, beta, data, X1, X2, kappa, W2_AA, MR_type, moment_type)
%{
 - Inputs - 
   delta:       nuisance parameter (to be estimated)
   beta:        strategic parameter 
   data:        data
   X1:          covariate for player 1
   X2:          covariate for player 2
   kappa:       tuning parameter for GMS
   W2_AA:       a vector of random variables used to implement the (multiplier) bootstrap.
   MR_type:     the type of resampling, i.e., DR or PR
   moment_type: whether aggregate or disaggregate moments

 - Outputs - 
   value:       function value 
%}

% decide if the number of moments are
n = size(data,1); % sample size

% the number of moment inequalities (which should appear first) and total moments;
if strcmp(moment_type, 'disagg')
    p = 8;
    k = 16;
elseif strcmp(moment_type, 'agg')
    p = 2;
    k = 4;
end

[~, mData, xi] = dataOp(delta, beta, data, kappa, X1, X2, moment_type);

if MR_type == 1 % DR method;
    value = S_function(W2_AA*zscore(mData, 1)/sqrt(n) + repmat(phi_function(xi',p,k-p),size(W2_AA,1),1), p);
elseif MR_type == 2 % PR method;
    value = S_function(W2_AA*zscore(mData, 1)/sqrt(n) + repmat(xi',size(W2_AA,1),1), p);
end
end
