% 5/30/2020 Shuowen Chen and Hiroaki Kaido
% Computes the Tn test statistics with Tikhonov regularizations
function [Test, g_delta] = stat(delta_hat, x, occurrence, n, lambda)
% Inputs:
%   delta_hat:  first-step estimated delta (1 by 2)
%   x:          possible combination of covariates (k by 2)
%   occurrence: number of occurrences of each combination of covariate and
%               potential outcomes (4k by 1)
%   n:          sample size
%   lambda:     regularization parameter
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
% need to adjust for this. 
% Two parts are regularized: varI and I_delta2 because the inverse of the
% latter is an input of varI

[~,D]=ldl(I_delta2); % LDL decomposition on I_delta2 
tol = 0.005;
if min(abs(diag(D))) < tol
    g_delta = Cbeta -  I_cross*(tikhonov(I_delta2,Cdelta,lambda));
    temp = I_cross';
    regsol = [tikhonov(I_delta2^(0.5),temp(:,1),lambda),...
        tikhonov(I_delta2^(0.5),temp(:,2),lambda)];
    varI = I_beta2 - regsol'*regsol;
else
    g_delta = Cbeta - I_cross/I_delta2 * Cdelta;
    varI = I_beta2 - I_cross/I_delta2 * I_cross';
end

[LI,DI] = ldl(varI);
DI(DI<tol) = tol; % truncate negative values

if min(abs(diag(DI))) <= tol
    Test = max(abs(tikhonov(LI*DI^(0.5),g_delta,lambda)));
else
    Test = max(abs(varI^(-0.5)*g_delta));
end
end
% Auxiliary function: Tikhonov regularization
function reg_sol = tikhonov(A,b,lambda)
    db = length(b);
    reg_sol = lsqr([A;lambda*eye(db)],[b;zeros(db,1)],[],1000);
end
