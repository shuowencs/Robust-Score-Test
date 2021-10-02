%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   5/1/2021 Shuowen Chen and Hiroaki Kaido
%   This is taken from code by Bugni, Canay and Shi (QE) with some changes
%   (1) We have two type of moment conditions
%   (2) p, the number of moment inequalities, are now in the fcn body
%   (3) we use dataOp.m to produce studendized moments
%   (4) instead of naming theta_to_min and theta_H0, we directly input delta
%       and beta

%   Defines sample criterion function according to MMM in Eq. (4.2) of BCS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = Qn_function_X(delta, beta, data, X1, X2, mtype)
%{
 - Inputs - 
   delta: nuisance parameter (to be estimated)
   beta:  strategic parameter 
   data:  data
   X1:    covariate for player 1
   X2:    covariate for player 2
   mtype: whether use aggregate or disaggregate moments

 - Outputs - 
   value: function value 
%}

% determines sample size;
n = size(data,1);

% the number of moment inequalities (which should appear first);
if strcmp(mtype, 'disagg')
    p = 8;
elseif strcmp(mtype, 'agg')
    p = 2;
end

% studentizes and averages the data (4 by S)
[mbar_std, ~, ~] = dataOp(delta, beta, data, 1, X1, X2, mtype); % No use for kappa, set to one.

% computes the sample criterion function;
value = S_function(sqrt(n)*mbar_std', p);
end
