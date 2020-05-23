% 5/20/2020 Shuowen Chen and Hiroaki Kaido
% Computes the Tn test statistics using the following
% regularization of the orthognalized variance covariance matrix
% varI = I_beta2 - kappa*I_cross/I_delta2 * I_cross';
function [Test, g_delta] = tnVar(delta_hat, x, occurrence, n, kappa)
% Inputs:
%   delta_hat:  first-step estimated delta (1 by 2)
%   x:          possible combination of covariates (k by 2)
%   occurrence: number of occurrences of each combination of covariate and
%               potential outcomes (4k by 1)
%   n:          sample size
%   kappa:      regularization parameter
% Outputs:
%   Test:       the sup test statistic
%   g_delta:    gn statistic
beta_null = [0, 0];
[zdelta, zbeta] = compute_z_all(beta_null, x, delta_hat);
Cbeta = zbeta * occurrence / sqrt(n);
Cdelta = zdelta * occurrence / sqrt(n);
% Compute the information matrix, which is approximated by the outer
% product of the scores. NOTE: we split the information matrix into 4
% parts, each of which is 2 by 2. 

% The following three blocks are each 4 by 16, where 16 stands for possible
% combinations of outcome s and covariate x. 
betablock = [zbeta(1,:).*zbeta(1,:); zbeta(1,:).*zbeta(2,:);...
    zbeta(2,:).*zbeta(1,:); zbeta(2,:).*zbeta(2,:)];

deltablock = [zdelta(1,:).*zdelta(1,:); zdelta(1,:).*zdelta(2,:);...
    zdelta(2,:).*zdelta(1,:); zdelta(2,:).*zdelta(2,:)];

crossblock = [zbeta(1,:).*zdelta(1,:); zbeta(1,:).*zdelta(2,:);...
    zbeta(2,:).*zdelta(1,:); zbeta(2,:).*zdelta(2,:)];

% The following three blocks are 4 by 1
betablock = betablock * occurrence / n;
deltablock = deltablock * occurrence / n;
crossblock = crossblock * occurrence / n;

% Reshape to be 2 by 2
I_beta2 = reshape(betablock, 2, 2);
I_cross = reshape(crossblock, 2, 2);
I_delta2 = reshape(deltablock, 2, 2);

% calculate the gn and sup test statistics
% Note: In small samples varI might contain negative variances in some 
% simulations. This causes trouble because the inverse of square root 
% contains complex numbers, and absolute value of complex numbers are real 
% numbers, which actually contaminates the test statistic. Therefore we 
% need to adjust for this. We consider regularizing the information 
% matrix to alleviate the issue. 
g_delta = Cbeta - I_cross/I_delta2 * Cdelta;
varI = I_beta2 - kappa*I_cross/I_delta2 * I_cross';
%varI = I_beta2 - I_cross/I_delta2 * I_cross';
%g_delta = Cbeta - I_cross/(I_delta2+kappa*[1,0;0,1]) * Cdelta;
%varI = I_beta2 - I_cross/(I_delta2+kappa*[1,0;0,1]) * I_cross';

% Test = max(abs(varI^(-0.5)*g_delta));
% a variable that indicates if varI^(-0.5) is real
temp = isreal(varI^(-0.5));
if temp == 1
    Test = max(abs(varI^(-0.5)*g_delta));
else
    Test = NaN;
end
end