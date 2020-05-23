% 5/21/2020 Shuowen Chen and Hiroaki Kaido
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
delta = [2, 2.5];
%% Standard deviations and finite-sample distribution of delta_hat
% We fix beta_alt and generate data of different sample sizes. We
% calculate the standard deviations of beta_hat in each sample size
% scenario (in each scenario we run MC S times).
N = [200, 500, 1000, 2000, 5000, 10000]; % five different sample sizes
S = 1000; % number of simulations for each sample size
% pre-allocation
std_N = zeros(2, length(N));
mean_N = zeros(2, length(N));
delta_hat_com_N = zeros(2*length(N), S);
for j = 1:length(N) % loop over sample sizes
    n = N(j);
    % Call data_generation function to generate data
    [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, ...
        delta, n, S, 'Complete');
    % Count Number of Occurrences of Events and X Configurations
    % There are in total 16 possible combinations of events and configurations.
    % Call function counting
    outcome_com = counting(data_com, x1_com, x2_com);
    % Conduct Restricted MLE
    delta_hat_com = zeros(2, S);
    for i = 1:S % loop over number of simulations
        % NOTE: rmle yields NA values while rmle2 doesn't 
        delta_hat_com(:,i) = rmle([0,0], [0,0], n, outcome_com(:,i));
        % delta_hat_com(:,i) = rmle2([0,0], [0,0], n, outcome_com(:,i));
    end
    % store the values
    delta_hat_com_N((2*j-1):(2*j),:) = delta_hat_com;
    % calculate the means and sample standard deviations
     mean_all = mean(delta_hat_com, 2, 'omitnan');
    %mean_all = mean(delta_hat_com, 2);
    % Row 1-2: delta1 and delta2 for iid; Row 3-4: delta1 and delta2 for non
    % iid; Row 5-6: delta1 and delta2 for LFP
    std_N(:,j) = sqrt(sum((delta_hat_com-mean_all).^2,2,'omitnan')./(n-1));
    %std_N(:,j) = sqrt(sum((delta_hat_com-mean_all).^2,2)./(n-1));
    mean_N(:,j) = mean_all;
end
% count number of NAs across sample sizes 
missing = sum(isnan(delta_hat_com_N),2); 
% Checking how close the distribution looks like normal
normplot(delta_hat_com_N(2,:)) % small sample is problematic, need to revise
normplot(delta_hat_com_N(11,:))
normplot(delta_hat_com_N(12,:)) % the tail is not satisfying
%% Size properties
% generate data under the null
size = 0.05; % nominal size
% critical value using the limiting distribution
cv = crt(size);
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
        delta_hat_com = rmle([0,0], [0,0], n, outcome_com(:,i));
        % generate var-cov matrices and gn statistics (generated under the null beta=0)
        [test(j,i), gn((2*j-1):(2*j),i)] = stat(delta_hat_com, X, outcome_com(:,i), n);
    end 
end
% count missing gn values, can see coincides with missings from delta_hat
gn1 = gn([1,3,5,7,9,11],:);
gn2 = gn([2,4,6,8,10,12],:);
missing2 = sum(isnan(gn1),2); 
% plot gn statistics
normplot(gn2(6,:)); % tail is a bit off
% compute the real size
[real, l] = actualsize(test, N, cv);
% From l2, we see that without regularization there are many missing value
% for smalle samples
%% The previous simulation uses delta_hat, now use true delta value
% generate data under the null
size = 0.05; % nominal size
% critical value using the limiting distribution
cv = crt(size);
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
% count missing gn values, can see coincides with missings from delta_hat
gn1_truedelta = gn_truedelta([1,3,5,7,9,11],:);
gn2_truedelta = gn_truedelta([2,4,6,8,10,12],:);
% plot gn_truedelta statistics
normplot(gn2_truedelta(6,:)); % tail looks good
% compute the real size
[real2, l2] = actualsize(test_truedelta, N, cv);
% From l2, we see that without regularization there are many missing value
% for smalle samples
%% Power analysis
% pre-allocation to store sup test statistic for each selection mechanism
test_iid = zeros(length(N), S);
test_non = zeros(length(N), S);
test_LFP = zeros(length(N), S);
k = 40; % local alternatives
for j = 1:length(N) % loop over sample size
    n = N(j);
    % Define true beta that generates the data (Change between beta_alt
    % and beta_alt2 in the data_generation function)
    beta_alt = h_alt(k,:)/sqrt(n);
    beta_alt2 = h_alt2(k,:)/sqrt(n);
    % Call data_generation function to generate data
    [data_iid, x1_iid, x2_iid] = data_generation(beta_alt, covariate, ...
        delta, n, S, selection1);
    [data_non, x1_non, x2_non] = data_generation(beta_alt, covariate, ...
        delta, n, S, selection2);
    [data_LFP, x1_LFP, x2_LFP] = data_generation(beta_alt, covariate, ...
        delta, n, S, selection3);
    outcome_iid = counting(data_iid, x1_iid, x2_iid);
    outcome_non = counting(data_non, x1_non, x2_non);
    outcome_LFP = counting(data_LFP, x1_LFP, x2_LFP);
    % loop over number of simulations
    for i = 1:S
        % Conduct Restricted MLE
        delta_hat_iid = rmle([0,0], [0,0], n, outcome_iid(:,i));
        delta_hat_non = rmle([0,0], [0,0], n, outcome_non(:,i));
        delta_hat_LFP = rmle([0,0], [0,0], n, outcome_LFP(:,i));
        % generate var-cov matrices and gn statistics (generated under the null beta=0)
        [test_iid(j,i), ~] = tn(delta_hat_iid, X, outcome_iid(:,i), n);
        [test_non(j,i), ~] = tn(delta_hat_non, X, outcome_non(:,i), n);
        [test_LFP(j,i), ~] = tn(delta_hat_LFP, X, outcome_LFP(:,i), n);
    end 
end
%% Now calculate the power of the test
power_iid = powerfcn(test_iid, N, crtvalue);
power_non = powerfcn(test_non, N, crtvalue);
power_LFP = powerfcn(test_LFP, N, crtvalue);
