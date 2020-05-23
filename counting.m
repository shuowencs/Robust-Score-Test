% 11/16/2019 Shuowen Chen and Hiroaki Kaido
% Counts the number of occurrences for each of the possible combination of 
% events and configurations of covariates. 
% Four potential outcomes: (0, 0), (1, 1), (1, 0) and (0, 1);
% Four covariate combinations: (1, 1), (1, -1), (-1, 1) and (-1, -1).
function outcome = counting(data, x1, x2)
% Inputs:
%   data: simulated observations (n by S) by the data_generation function
%   x1:   simulated x1 (n by S) by the data_generation function
%   x2:   simulated x2 (n by S) by the data_generation function
% Outputs:
%   outcome: a 16 by S matrix that lists the number of occurrences for each
%            Monte Carlo simulation sequence
%   note: S denotes number of simulations for each sample sizetemp = size(data);
temp = size(data);
S = temp(2);
% placeholder
outcome = NaN(16, S);
outcome(1, :) = sum(data(:, :)==0 & x1(:, :)==1 & x2(:, :)==1);
outcome(2, :) = sum(data(:, :)==11 & x1(:, :)==1 & x2(:, :)==1);
outcome(3, :) = sum(data(:, :)==10 & x1(:, :)==1 & x2(:, :)==1);
outcome(4, :) = sum(data(:, :)==1 & x1(:, :)==1 & x2(:, :)==1);
outcome(5, :) = sum(data(:, :)==0 & x1(:, :)==1 & x2(:, :)==-1);
outcome(6, :) = sum(data(:, :)==11 & x1(:, :)==1 & x2(:, :)==-1);
outcome(7, :) = sum(data(:, :)==10 & x1(:, :)==1 & x2(:, :)==-1);
outcome(8, :) = sum(data(:, :)==1 & x1(:, :)==1 & x2(:, :)==-1);
outcome(9, :) = sum(data(:, :)==0 & x1(:, :)==-1 & x2(:, :)==1);
outcome(10, :) = sum(data(:, :)==11 & x1(:, :)==-1 & x2(:, :)==1);
outcome(11, :) = sum(data(:, :)==10 & x1(:, :)==-1 & x2(:, :)==1);
outcome(12, :) = sum(data(:, :)==1 & x1(:, :)==-1 & x2(:, :)==1);
outcome(13, :) = sum(data(:, :)==0 & x1(:, :)==-1 & x2(:, :)==-1);
outcome(14, :) = sum(data(:, :)==11 & x1(:, :)==-1 & x2(:, :)==-1);
outcome(15, :) = sum(data(:, :)==10 & x1(:, :)==-1 & x2(:, :)==-1);
outcome(16, :) = sum(data(:, :)==1 & x1(:, :)==-1 & x2(:, :)==-1);
end