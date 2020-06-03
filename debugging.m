% 5/29/2020 Shuowen Chen and Hiroaki Kaido
% This file checks the properties of 
% (1) Estimated nuisance parameters
%   Problem:
%     For small samples, current algorithm produces NA values in some simulations. 
%     In these cases, I alternated to steepest ascent and managed to get estimands,
%     but the empirical distribution of estimated parameters looks strange.
%     
%   Current conjecture:
%     needs some refinement in the algorithms (stepsize), essentially the
%     same problem as in calculating the sup test statistic
% 
% (2) The sup test statistic
%   Problem:
%     The test statstic produces NA or complex values in some simulations,
%     even for modeate number of sample sizes. Two sources of problems:
%     (1) NA delta estimates carries over to produce NA test statstic
%     (2) The gn statistic contains the inverse of I_deltadelta matrix
%     the sup test statistic contains the inverse of square root of the
%     following matrix:
%     I_betabeta - I_betadelta * inv(I_deltadelta) * I_deltabeta
%     In finite sample can lead to complex values
%   Current conjecture: 
%     simulations that generate extremely large test statstic values stem from
%     bad scaling of the inverse matrix. Need to substantiate this conjecture
%     and fix it: more careful choice of regularization parameters 
%   We first check the performance under the null (complete model)
%% Setting Parameters            
% beta null
beta0 = 0;         
% In each sample we will draw X from this vector uniformly with replacement
covariate = [1, -1];
% true value of nuisance parameters
delta = [2, 1.5];
%% Standard deviations and finite-sample distribution of delta_hat
% We fix beta_alt and generate data of different sample sizes. We
% calculate the standard deviations of beta_hat in each sample size
% scenario (in each scenario we run MC S times).
% We prefer the LM-BFGS algorithm. The code is provided by Stephen Becker
N = [200, 500, 1000, 2000, 5000, 10000]; % six different sample sizes
S = 5000; % number of simulations for each sample size
% pre-allocation
%delta_BHHH = zeros(2*length(N), S);
%delta_BFGS = zeros(2*length(N), S);
delta_LM = zeros(2*length(N), S);
%est1 = zeros(2, S);
%est2 = zeros(2, S);
est3 = zeros(2, S);
%parpool(8) % for parallelization (adjust this depending on # of cores)
for j = 1:length(N) % loop over sample sizes
    n = N(j);
    % Call data_generation function to generate data
    [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, ...
        delta, n, S, 'Complete');
    % loop over simulations
    for i = 1:S 
        % Conduct Restricted MLE
        %est1(:,i) = rmle([0,0], n, 'BHHH', data_com(:,i), x1_com(:,i), x2_com(:,i));
        %est2(:,i) = rmle([0,0], n, 'BFGS', data_com(:,i), x1_com(:,i), x2_com(:,i));
        est3(:,i) = rmle([0,0], n, 'LM-BFGS', data_com(:,i), x1_com(:,i), x2_com(:,i));
    end
    % Store estimates
    %delta_BHHH((2*j-1):(2*j), :) = est1;
    %delta_BFGS((2*j-1):(2*j), :) = est2;
    delta_LM((2*j-1):(2*j), :) = est3;
end
%% Checking property for results (by BHHH)
% count number of NAs across sample sizes 
missing_BHHH = sum(isnan(delta_BHHH),2); 
missing_BFGS = sum(isnan(delta_BFGS),2);
missing_LM = sum(isnan(delta_LM),2);
% Checking how close the distribution looks like normal
normplot(delta_BHHH(2,:)) % small sample is problematic, need to revise
normplot(delta_BHHH(11,:))
normplot(delta_BHHH(12,:)) % the tail is not satisfying
J = length(N);

figure(1)
for j=1:J
    subplot(J,2,2*j-1)
    histogram(delta_BHHH(2*(j-1)+1,:),floor(2*S^(1/3)));
    title(['n=' num2str(N(j))])
    subplot(J,2,2*j)
    histogram(delta_BHHH(2*j,:),floor(2*S^(1/3)));
    title(['n=' num2str(N(j))])
end
%% Size properties
% generate data under the null
level = 0.05; % nominal size
% critical value using the limiting distribution
cv = crt(level);
N = [200, 500, 1000, 1500, 20000, 30000];
S = 1000; % number of simulations in each DGP
X = [1,1; 1,-1; -1,1; -1,-1]; % combination of covariates
% pre-allocation to store g_delta statstics and sup test statistic
gn = zeros(2*length(N), S);
test = zeros(length(N), S);
for j = 1:length(N) % loop over sample size
    n = N(j);
    % Call data_generation function to generate data
    [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, ...
        delta, n, S, 'Complete');
    outcome_com = counting(data_com, x1_com, x2_com);
    % loop over number of simulations
    for i = 1:S
        % Conduct Restricted MLE
        delta_hat_com = rmle([0,0], n, 'LM-BFGS', data_com(:,i), x1_com(:,i), x2_com(:,i));
        % generate var-cov matrices and gn statistics (generated under the null beta=0)
        [test(j,i), gn((2*j-1):(2*j),i)] = stat(delta_hat_com, X, outcome_com(:,i), n, 0.1);
    end 
end
% compute the real size
[real, l] = actualsize(test, N, cv);
%% The previous simulation uses delta_hat, now use true delta value
% generate data under the null
level = 0.05; % nominal size
% critical value using the limiting distribution
cv = crt(level);
N = [200, 500, 1000, 1500, 20000, 30000];
S = 1000; % number of simulations in each DGP
X = [1,1; 1,-1; -1,1; -1,-1]; % combination of covariates
% pre-allocation to store g_delta statstics and sup test statistic
gn_truedelta = zeros(2*length(N), S);
test_truedelta = zeros(length(N), S);

for j = 1:length(N) % loop over sample size
    n = N(j);
    % Call data_generation function to generate data
    [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, ...
        delta, n, S, 'Complete');
    outcome_com = counting(data_com, x1_com, x2_com);
    % loop over number of simulations
    for i = 1:S
        % using the true delta
        [test_truedelta(j,i), gn_truedelta((2*j-1):(2*j),i)] = stat(delta, X, outcome_com(:,i), n);
    end 
end
% compute the real size
[real2, l2] = actualsize(test_truedelta, N, cv);
%% Power analysis
N = [200, 500, 1000];
% pre-allocation to store sup test statistic for each selection mechanism
test_iid = zeros(length(N), S);
test_non = zeros(length(N), S);
test_LFP = zeros(length(N), S);
% pre-define alternatives
h_alt  = -(0:0.1:5)';
h_alt2 = -(0:0.1:5)'; 
k = 40; % local alternatives
for j = 1:length(N) % loop over sample size
    n = N(j);
    % Define true beta that generates the data (Change between beta_alt
    % and beta_alt2 in the data_generation function)
    beta_alt = h_alt(k,:)/sqrt(n);
    beta_alt2 = h_alt2(k,:)/sqrt(n);
    % Call data_generation function to generate data
    [data_iid, x1_iid, x2_iid] = data_generation([beta_alt,beta_alt2], covariate, ...
        delta, n, S, 'IID');
    [data_non, x1_non, x2_non] = data_generation([beta_alt,beta_alt2], covariate, ...
        delta, n, S, 'Non IID');
    [data_LFP, x1_LFP, x2_LFP] = data_generation([beta_alt,beta_alt2], covariate, ...
        delta, n, S, 'LFP');
    outcome_iid = counting(data_iid, x1_iid, x2_iid);
    outcome_non = counting(data_non, x1_non, x2_non);
    outcome_LFP = counting(data_LFP, x1_LFP, x2_LFP);
    % loop over number of simulations
    for i = 1:S
        % Conduct Restricted MLE
        delta_hat_iid = rmle([0,0], n, 'LM-BFGS', data_iid(:,i), x1_iid(:,i), x2_iid(:,i));
        delta_hat_non = rmle([0,0], n, 'LM-BFGS', data_non(:,i), x1_non(:,i), x2_non(:,i));
        delta_hat_LFP = rmle([0,0], n, 'LM-BFGS', data_LFP(:,i), x1_LFP(:,i), x2_LFP(:,i));
        % generate var-cov matrices and gn statistics (generated under the null beta=0)
        [test_iid(j,i), ~] = stat(delta_hat_iid, X, outcome_iid(:,i), n, 0.1);
        [test_non(j,i), ~] = stat(delta_hat_non, X, outcome_non(:,i), n, 0.1);
        [test_LFP(j,i), ~] = stat(delta_hat_LFP, X, outcome_LFP(:,i), n, 0.1);
    end 
end
%% Now calculate the power of the test
power_iid = powerfcn(test_iid, N, cv);
power_non = powerfcn(test_non, N, cv);
power_LFP = powerfcn(test_LFP, N, cv);