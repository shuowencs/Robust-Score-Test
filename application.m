% 10/1/2021 Shuowen Chen and Hiroaki Kaido
% Application to data in Kline and Tamer (2016, QE)
% We do not consider correlation between errors
%% Data Preparation (With Discretization)
data = dlmread('airlinedata.dat');
% sample size
n = size(data, 1);
% Outcome variable: LCC first and OA second
y = 10*data(:, 4) + data(:, 5); 
% Construct regressors. Each player has three regressors:
% own presence (col1), market size (col2) and constant (col3)
x_LCC = ones(n, 3);
x_OA = ones(n, 3);
% Two excluded regressors
LCCmedianPres = median(data(:, 1));
OAmedianPres = median(data(:, 2));
x_LCC(:, 1) = (data(:, 1) >= LCCmedianPres);
x_OA(:, 1) = (data(:, 2) >= OAmedianPres);
% One common regressors
sizeMedian = median(data(:, 3));
x_LCC(:, 2) = (data(:, 3) >= sizeMedian);
x_OA(:, 2) = (data(:, 3) >= sizeMedian);
%% Setting Parameter Space for Searching and Initial Guesses
% Parameter space of each parameter (based on figure 2 in Kline and Tamer, 2016)
% order of parameter: LCC presence, size and constant, OA presence, size
% and constant
lb = [-0.5; -0.25; -2.5; -0.25; -0.1; -1.5];
ub = [3.5; 1; 1.5; 1.5; 1.25; 1];
% Generate multiple staring points
M = 150;
initial = lhsscale(M, 6, ub, lb);
%% RMLE (With Discretization)
% For each starting point, compute delta_est
delta_pool = zeros(6, M);
fvalue = zeros(M, 1);
for i = 1:M
    [delta_pool(:,i), fvalue(i)] = rmleapp(initial(i,:), 'BFGS', y, x_LCC, x_OA, lb, ub);
end
% Find the parameters that yield the smallest objective function (largest log likelihood)
indice = find(fvalue == min(fvalue));
delta_est = delta_pool(:, indice);
%% Test statistics (With Discretization)
[test_stat, g_delta, varI, h, Vstar] = statapp(delta_est, y, x_LCC, x_OA);
% Simulate critical values (this step takes some time)
crit05 = crt_ab(0.05, varI);
crit01 = crt_ab(0.01, varI);
<<<<<<< HEAD
% Display test result
disp(['=== Test Statistic is ', num2str(test_stat), ...
    ', 95th percentile critical value is ', num2str(crit05), ...
    ', 99th percentile critical value is ', num2str(crit01)]); 
% Compare with Shrinkage CV
=======
pval   = pval_ab(test_stat,varI);

%% Display test result
disp(['=== Test Statistic is ', num2str(test_stat), ...
    ', 95th percentile critical value is ', num2str(crit05), ...
    ', 99th percentile critical value is ', num2str(crit01),...
    ', p-value is ', num2str(pval)]); 
%% Implement post model inference
>>>>>>> 9eef4acfa702b06afb13b459c913acc33c82e1f4
% (first add KMS folder to the path)
kappa = 1/sqrt(log(n)); % shrinkage factor
disp(['=== Test Statistic is ', num2str(test_stat), ...
    ', kappa*cv05 is ', num2str(kappa*crit05), ...
    ', kappa*cv01 is ', num2str(kappa*crit01), ...
    '. Therefore KMS should be implemented.']);
%% Redux without Discretization
% Construct regressors. Each player has three regressors:
% own presence (col1), market size (col2) and constant (col3)
x2_LCC = data(:, [1, 3]);
x2_OA = data(:, [2, 3]);
% We need to rescale the covariate because the common regressor is huge
colmax = max(x2_LCC);
colmin = min(x2_LCC);
x2_LCC = rescale(x2_LCC, 'InputMin',colmin,'InputMax',colmax);
colmax = max(x2_OA);
colmin = min(x2_OA);
x2_OA = rescale(x2_OA, 'InputMin',colmin,'InputMax',colmax);
% add constant
x2_LCC = [x2_LCC, ones(n, 1)];
x2_OA = [x2_OA, ones(n, 1)];
%% RMLE (Without Discretization)
lb = [-0.5; -0.25; -2.5; -0.25; -0.1; -1.5];
ub = [3.5; 1; 1.5; 1.5; 1.25; 1];
% Generate multiple staring points
M = 150;
initial = lhsscale(M, 6, ub, lb);
delta_pool_nodis = zeros(6, M);
fvalue2 = zeros(M, 1);
for i = 1:M
    [delta_pool_nodis(:,i), fvalue2(i)] = rmle_nodis(initial(i,:), 'LM-BFGS', y, x2_LCC, x2_OA, lb, ub);
end
% Find the parameters that yield the smallest objective function (largest log likelihood)
indice2 = find(fvalue2 == min(fvalue2));
delta_est_nodis = delta_pool_nodis(:, indice2);
%% Simulate critical values (this step takes some time)
[test_nodis, ~, var_nodis, ~, ~] = test_stat_no_discretization(delta_est_nodis, y, x2_LCC, x2_OA);
% Simulate critical values (this step takes some time)
crit05_nodis = crt_ab(0.05, var_nodis);
crit01_nodis = crt_ab(0.01, var_nodis);
%% Display test result
disp(['=== Test Statistic without discretization is ', num2str(test_nodis), ...
    ', 95th percentile critical value is ', num2str(crit05_nodis), ...
    ', 99th percentile critical value is ', num2str(crit01_nodis)]); 
kappa = 1/sqrt(log(n)); % shrinkage factor
disp([', kappa*cv05 is ', num2str(kappa*crit05_nodis), ...
    ', kappa*cv01 is ', num2str(kappa*crit01_nodis), ...
    '. Therefore KMS should be implemented.']);
