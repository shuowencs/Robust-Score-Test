% 6/16/2020 Shuowen Chen and Hiroaki Kaido
% Conducts the restricted mle
% Under the restriction beta1 = beta2 = 0.
function [delta] = rmle(delta_initial, n, algorithm, data, x1, x2)
% Input:
%   delta_initial:  initial guess of coefficients of covariates (1 by 2*d)
%   n:              sample size
%   algorithm:      the algorithm for estimation: BFGS, LM-BFGS or BHHH
%   data:           game outcome
%   x1:             covariates of player 1, same definition for x2 (n by d)
% Output: 
%   delta:          estimates of parameters (1 by 2*d)

% Compute the outcome for each configuration combination of x and delta
% The 16 by S matrix that stores the number of occurrences of each 
% combination of events and covariate configurations
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
    f = @(delta_est) obj(delta_est, outcome);
    % Setting Options
    % Once gradient is supplied, switch to trust region algorithm
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
    %initial = [1;1];
    % Call fminunc
    delta = fminunc(f, delta_initial', options);
elseif strcmp(algorithm, 'LM-BFGS')
    % Now use the LM-BFGS algorithm
    fcn = @(nuisance) loglikelihood(nuisance, outcome);
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
% Computes the least favorable pairs under the null
% -- Input --  
%  delta:   nuisance parameter (for estimation purpose): 2 by 1 
% -- Output --
%  Q0:     Distribution under the null [q0_00, q0_11, q0_10, q0_01]
function Q0 = glfp(delta)
x = [1,1;1,-1;-1,1;-1,-1];
xdelta1 = x(:,1)*delta(1);
xdelta2 = x(:,2)*delta(2);
Phi1 = normcdf(xdelta1);
Phi2 = normcdf(xdelta2);
% Each is 4 by 1
q0_00 = (1-Phi1).*(1-Phi2);
q0_11 = Phi1.*Phi2;
q0_10 = Phi1.*(1-Phi2);
q0_01 = (1-Phi1).*Phi2;
Q0 = [q0_00', q0_11', q0_10', q0_01'];
end

% Calculate the nagtive loglikelihood and corresponding score for BFGS
% outcome: frequencies of covariate and outcome combinations (16 by 1)
function [f, g] = obj(delta, outcome)
q0 = glfp(delta');
frequency = reshape(outcome, 4, 4);
frequency = reshape(frequency', 1, 16);
f = -sum(sum(frequency .* log(q0)));
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

%%%% For LM-BFGS
% Calculate the negative of loglikelihood function
function f = loglikelihood(delta, outcome)
% delta should be 2 by 1
q0 = glfp(delta');
frequency = reshape(outcome, 4, 4);
frequency = reshape(frequency', 1, 16);
f = -sum(sum(frequency .* log(q0)));
end

% Calculate the score of the negative loglikelihood function
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