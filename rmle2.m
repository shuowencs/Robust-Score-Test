% 4/20/2020 Shuowen Chen and Hiroaki Kaido
% A temp file that refines the algorithm
% Conducts the restricted mle via BHHH algorithm
% Under the restriction beta1 = beta2 = 0.
function [delta] = rmle2(beta, delta_initial, n, outcome)
% Input:
%   beta:           structural parameters of interest under restriction (0,0)
%   delta_initial:  initial guess of coefficients of covariates (1 by 2)
%   n:              sample size
%   outcome:        the 16 by S matrix that stores the number of
%                   occurrences of each combination of events and covariate
%                   configurations
% Output: 
%   delta:          estimates of parameters (1 by 2)

% Acknowledgement: 
% tolerance level and statistics for convergence evaluation are taken from
% Chapter 8, Train, Kenneth, Discrete Choice Methods with Simulation, 2009,
% 2nd edition, Cambridge University Press. 
tol = 0.0001;
stat = 1;
% stepsize for BHHH algorithm
lambda = 1;

while stat > tol
    % Call compute_z_all to get zdelta, 2 by 16. 
    [zdelta, ~] = compute_z_all(beta, [1,1; 1,-1; -1,1; -1,-1], delta_initial);
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
    delta_update = delta_initial' + lambda*hess\sum_l_dot;
    delta_initial = delta_update';
    % Update statistic m to evaluate convergence
    stat = sum_l_dot'*(hess\sum_l_dot)/n;
end
delta = delta_initial;

% For some simulations in small samples, the BHHH algorithm provides NaN
% values as the inverse of the outer product of the scores becomes singular
% or badly scaled. For these cases, we use the steepest ascent method as an
% alternative
if sum(isnan(delta)) == 2
    lambda2 = 0.001; % stepsize 
    tol = 0.0001;
    stat = 1;
    delta_initial = [0, 0];
    I2 = [1,0;0,1];
    while stat > tol
        % Call compute_z_all to get zdelta, 2 by 16.
        [zdelta, ~] = compute_z_all(beta, [1,1; 1,-1; -1,1; -1,-1], delta_initial);
        % Compute sum of score functions (2 by 1)
        sum_l_dot = zdelta * outcome;
        % Steepest ascent update (NOTE delta_update is 2 by 1)
        delta_update = delta_initial' + lambda2*I2*sum_l_dot;
        delta_initial = delta_update';
        % Update statistic m to evaluate convergence
        stat = sum_l_dot'*I2*sum_l_dot/n;
    end
    delta = delta_initial;
end
end