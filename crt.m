% 3/10/2020 Hiroaki Kaido and Shuowen Chen
% The function calculates the critical value of the tn statistic
% Recall that the limiting distribution of tn is the sup of bivariate
% standard normal distribution
function crit = crt(nominal)
    rng(123); % set seed for random number generations
    J = 2;    % # of components in the score
    R = 5000; % # of draws
    z = randn(J, R); % J by R vector of draws from standard normal
    sim_stat = max(abs(z)); % limiting distribution
    crit = quantile(sim_stat, (1-nominal));
end