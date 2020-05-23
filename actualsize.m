% 5/20/2020
% Shuowen Chen and Hiroaki Kaido
% A function that calculates the critical values of the simulated Tn
% statistics
function [real, l] = actualsize(test_matrix, N, cv)
% Outputs:
% real: the actual size
% l: (auxiliary output) length of sequence (some simulations end up producing NA values, 
%    so compute l to have a sense how severe the issue is)
real = zeros(length(N), 1); % actual size placeholder
l = zeros(length(N), 1);
for i = 1:length(N)
    test_seq = test_matrix(i, :);
    % this is to delete ill-behaved test statistics (due to var-cov matrix)
    % test_seq = test_seq(~isnan(test_seq) & imag(test_seq) == 0);
    test_seq = test_seq(~isnan(test_seq));
    % compute the actual size
    real(i) = sum(test_seq > cv)/length(test_seq);
    % compute the length of sequence
    l(i) = length(test_seq);
end
end