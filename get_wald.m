
function [cil,ciu] =  get_wald(delta_hat,x,occurrence,n,level,comp)
% Inputs
% delta_hat: 1st stage estimator
% x:          support of covariates
% occurrence: outcome frequency
% n:          sample size
% level:      nominal level
% comp:       component of interest (=1)
beta_null = [0, 0];
[zdelta, ~] = compute_z_all(beta_null, x, delta_hat);
deltablock = [zdelta(1,:).*zdelta(1,:); zdelta(1,:).*zdelta(2,:);...
    zdelta(2,:).*zdelta(1,:); zdelta(2,:).*zdelta(2,:)];
deltablock = deltablock * occurrence / n;
I_delta2 = reshape(deltablock, 2, 2);
I_delta2 = get_sigmatilde(I_delta2); % regularize using Andrews & Barwic
Sigma = inv(I_delta2);
c = norminv(1-level/2,0,1);
cil = delta_hat(comp)-c*sqrt(Sigma(comp,comp)/n);
ciu = delta_hat(comp)+c*sqrt(Sigma(comp,comp)/n);
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