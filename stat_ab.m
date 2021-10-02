% 5/21/2020 Shuowen Chen and Hiroaki Kaido
% Computes the Tn test statistics without regularizations
function [Test, g_delta, varI] = stat_ab(delta_hat, x, occurrence, n)
% Inputs:
%   delta_hat:  first-step estimated delta (1 by 2)
%   x:          possible combination of covariates (k by 2)
%   occurrence: number of occurrences of each combination of covariate and
%               potential outcomes (4k by 1)
%   n:          sample size
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

I_delta2 = get_sigmatilde(I_delta2); % regularize using Andrews & Barwick
g_delta = Cbeta - I_cross/I_delta2 * Cdelta;
varI = I_beta2 - I_cross/I_delta2 * I_cross';
varI = get_sigmatilde(varI);

gtilde = varI^(-0.5)*g_delta;
Test = gtilde'*gtilde;
Vmin = @(x) quadform(varI,g_delta-x);
[~,Vstar] = fmincon(Vmin,[0;0],eye(2),zeros(2,1));
Test = Test - Vstar;
end

function SigmaTilde = get_sigmatilde(Sigma)
% Inputs:
% Sigma: Covariance matrix
% SigmaTilde: Regularized covariance matrix using Andrews & Barwick's
% method
sh_param = 0.05;
sigma = sqrt(diag(Sigma));
Omega  = Sigma./(sigma*sigma');
D = diag(diag(Sigma));
SigmaTilde = Sigma + max(sh_param-det(Omega),0).*D;
end

function t = quadform(A,g)
gtilde = A^(-0.5)*g;
t = gtilde'*gtilde;
end