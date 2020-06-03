% 5/28/2020 Shuowen Chen and Hiroaki Kaido
% Conducts the restricted mle
% Under the restriction beta1 = beta2 = 0.
function [delta] = rmle(delta_initial, n, algorithm, data, x1, x2)
% Input:
%   beta:           structural parameters of interest under restriction (0,0)
%   delta_initial:  initial guess of coefficients of covariates (1 by 2)
%   n:              sample size
%   outcome:        the 16 by S matrix that stores the number of
%                   occurrences of each combination of events and covariate
%                   configurations
%   algorithm:      the algorithm for estimation: BFGS, LM-BFGS or BHHH
% Output: 
%   delta:          estimates of parameters (1 by 2)

% Compute the outcome for each configuration combination of x and delta
outcome = counting(data, x1, x2);
if strcmp(algorithm, 'BHHH')
    % Acknowledgement:
    % tolerance level and statistics for convergence evaluation are taken from
    % Chapter 8, Train, Kenneth, Discrete Choice Methods with Simulation, 2009,
    % 2nd edition, Cambridge University Press.
    tol = 0.0001;
    stat = 1;
    % stepsize for BHHH algorithm
    lambda = 0.8;
    % beta under the restriction
    beta = [0,0];
    while stat > tol
        % Call compute_z_all to get zdelta, 2 by 16.
        [zdelta, ~] = compute_z_all(beta, [1,1; 1,-1; -1,1; -1,-1], delta_initial');
        zdelta_t = zdelta';
        % Compute sum of score functions (2 by 1)
        sum_l_dot = zdelta * outcome;
        % Compute outer product of score functions
        hess = zeros(2, 2);
        for i = 1:16
            score_sqr_i = outcome(i)*zdelta(:, i)*zdelta_t(i, :);
            hess = hess + score_sqr_i;
        end
        % BHHH update (NOTE delta_update is 2 by 1)
        delta_update = delta_initial' + lambda*(hess+log(n)/sqrt(n)*[1,0;0,1])\sum_l_dot;
        %delta_update = delta_initial' + lambda*hess\sum_l_dot;
        delta_initial = delta_update';
        % Update statistic m to evaluate convergence
        stat = sum_l_dot'*((hess+log(n)/sqrt(n)*[1,0;0,1])\sum_l_dot)/n;
        %stat = sum_l_dot'*(hess\sum_l_dot)/n;
    end
    delta = delta_initial;
elseif strcmp(algorithm, 'BFGS')
    % Objective function and gradient
    f = @(delta_est) obj(data, x1, x2, delta_est, outcome);
    % Setting Options
    % Once gradient is supplied, switch to trust region algorithm
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
    %initial = [1;1];
    % Call fminunc
    delta = fminunc(f, delta_initial', options);
elseif strcmp(algorithm, 'LM-BFGS')
    % Now use the LM-BFGS algorithm
    fcn = @(nuisance) loglikelihood(data, x1, x2, nuisance);
    % Gradient
    grad = @(nuisance) score(nuisance, outcome);
    l  = zeros(2,1);    % lower bound
    u  = inf(2,1);         % upper bound
    fun = @(nuisance) fminunc_wrapper(nuisance, fcn, grad);
    % Request very high accuracy for this test:
    opts2 = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
    %opts2.printEvery = 100;
    % Run the algorithm:
    [delta, ~, ~] = lbfgsb(fun, l, u, opts2);
end
end

% --- Auxiliary Functions --- 

% 5/26/2020 Shuowen Chen and Hiroaki Kaido
% Computes the neg log likelihood function and gradient
% -- Inputs --  
%  data:    market outcome (n by 1)
%  x1:      covariate for player 1 (n by 1)
%  x2:      covariate for player 2 (n by 1)
%  delta:   nuisance parameter (for estimation purpose): 2 by 1 
%  outcome: frequencies of covariate and outcome combinations (16 by 1)
% -- Outputs --
%  f:       negative log likelihood function
%  g:       score function of f
function [f, g] = obj(data, x1, x2, delta, outcome)
% compute the PMF for each market given covariates and delta (n by 4):
% [q0_00, q0_11, q0_10, q0_01]
q0 = get_LFP(x1, x2, delta');
% construct the negative log likelihood function
temp1 = zeros(size(data));
temp2 = zeros(size(data));
temp3 = zeros(size(data));
temp4 = zeros(size(data));
temp1(data == 0) = 1;
temp2(data == 11) = 1;
temp3(data == 10) = 1;
temp4(data == 1) = 1;
temp = [temp1, temp2, temp3, temp4];
% objective function
f = -sum(sum(temp .* q0));
% gradient required
if nargout > 1 
    % Call compute_z_all to get zdelta, 2 by 16.
    [zdelta, ~] = compute_z_all([0, 0], [1,1; 1,-1; -1,1; -1,-1], delta);
    % Compute sum of score functions (2 by 1)
    g = zdelta * outcome;
    % Since we add negative sign to use minimization algorithm, adjust the
    % score accordingly
    g = -g;
end
end

function Q0 = get_LFP(x1, x2, delta)
% Purpose: The function calculates the analytical LFP given and Xdelta. 
% Inputs:
%   beta:  structural parameters of interest
%   x1:    covariate of player 1; n by 1 matrix
%   x2:    covariate of player 2; n by 1 matrix
%   delta: nuisance parameter; 1 by 2
% Outputs:
%   Q0:     Distribution under the null [q0_00, q0_11, q0_10, q0_01],
%           Dimension: nS by 4. 
xdelta1 = x1*delta(1);
xdelta2 = x2*delta(2);
dim = size(xdelta1);
% sample size
n = dim(1);
% number of MC
S = dim(2);
% Each of the following four objects are n by S. 
Phi1 = normcdf(xdelta1);
Phi2 = normcdf(xdelta2);
% Within a MC, Q0 is identical arocss all subcases. 
% The Q0 below is a n*S by 4 matrix 
Q0 = [reshape((1-Phi1).*(1-Phi2),n*S,1), reshape(Phi1.*Phi2,n*S,1),...
    reshape(Phi1.*(1-Phi2),n*S,1), reshape((1-Phi1).*Phi2,n*S,1)];
end

% for LM-BFGS
function f = loglikelihood(data, x1, x2, delta)
% delta should be 2 by 1
q0 = get_LFP(x1, x2, delta');
% construct the negative log likelihood function
temp1 = zeros(size(data));
temp2 = zeros(size(data));
temp3 = zeros(size(data));
temp4 = zeros(size(data));
temp1(data == 0) = 1;
temp2(data == 11) = 1;
temp3(data == 10) = 1;
temp4(data == 1) = 1;
temp = [temp1, temp2, temp3, temp4];
% objective function
f = -sum(sum(temp .* q0));
end

function g = score(delta, outcome)
% delta should be 2 by 1
% Call compute_z_all to get zdelta, 2 by 16.
[zdelta, ~] = compute_z_all([0, 0], [1,1; 1,-1; -1,1; -1,-1], delta);
% Compute sum of score functions (2 by 1)
g = zdelta * outcome;
% Since we add negative sign to use minimization algorithm, adjust the
% score accordingly
g = -g;
end