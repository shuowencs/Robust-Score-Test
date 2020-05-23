% 5/21/2020
% Shuowen Chen and Hiroaki Kaido
% A function that computes the power of the Tn test statistic
function [power, l] = powerfcn(test_matrix, N, cv)
% inputs: 
% test_matrix: matrix that stores the simulated test statistics across N
% N: vector of sample sizes
% cv: critical value based on limiting distribution
power = zeros(length(N), 1); % placeholder
for i = 1:length(N)
    test_seq = test_matrix(i, :);
    l = length(test_seq);
    test_seq = test_seq(~isnan(test_seq));
    power(i) = 1 - sum(test_seq < cv)/length(test_seq);
end
end