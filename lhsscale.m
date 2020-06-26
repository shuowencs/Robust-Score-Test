% 6/17/2020 Shuowen Chen and Hiroaki Kaido
% Implements Latin Hypercube Sampling and rescale to produce initial 
% guesses for restricted MLE
function output = lhsscale(nsamples, nvars, ux, lx)
% Inputs:
%   nsamples: the number of points for each parameter
%   nvars:    the dimension of the parameter
%   ux:       the upper bound of the parameter space (1 by nvars)
%   lx:       the lower bounds of the parameter space (1 by nvars)
rng(123)
% draw points on a cube (Latin hypercube sampling)
% the nsample points will be from (0,1/n), (1/n, 2/n),...,(1-1/n,1), where
% n is shorthand notation of nsamples. Note the intervals can be randomly
% permutated
temp = lhsdesign(nsamples, nvars);  
% rescale to draw points on parameter space	
output = zeros(nsamples, nvars); 
for i = 1:nvars
    output(:, i) = repmat((ux(i)-lx(i)), nsamples, 1).*temp(:,i) + ...
        repmat(lx(i), nsamples, 1);
end
end

