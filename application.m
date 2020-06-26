% 6/16/2020 Shuowen Chen and Hiroaki Kaido
% Application of robust score test to data in Kline and Tamer (2016, QE)
%% Data Preparation
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
%% We first shut down correlation between errors
% Parameter space of each parameter (based on figure 2 in Kline and Tamer, 2016)
% order of parameter: LCC presence, size and constant, OA presence, size
% and constant
lb = [-0.5; -0.25; -2.5; -0.25; -0.1; -1.5];
ub = [3.5; 1; 1.5; 1.5; 1.25; 1];
% Generate multiple staring points
M = 50;
initial = lhsscale(M, 6, ub, lb);
%% RMLE
% For each starting point, compute delta_est
delta_pool = zeros(6, M);
fvalue = zeros(M, 1);
for i = 1:M
    [delta_pool(:,i), fvalue(i)] = rmleapp(initial(i,:), 'LM-BFGS', y, x_LCC, x_OA, lb, ub);
end
% Find the parameters that yield the smallest objective function (largest log likelihood)
indice = find(fvalue == min(fvalue));
delta_est = delta_pool(:, indice);
%% Sup test statistics
test_stat = statapp(delta_est, y, x_LCC, x_OA);
crit05 = crt(0.05);
crit01 = crt(0.01);