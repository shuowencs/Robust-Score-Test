% 10/1/2021 Shuowen Chen and Hiroaki Kaido
% Implements restricted MLE for application without discretization
function [delta, fval] = rmle_nodis(delta_initial, algorithm, data, x1, x2, lb, ub)
% Input:
%   delta_initial: Initial guess of delta (1 by 6)
%   algorithm:     Algorithm for estimation 
%   data:          market outcomes   
%   x1:            covariates of player 1, same definition for x2 (n by d)
%   lb:            If algorithm is 'LM-BFGS', specify the lower bound of
%                  parameter space
%   ub:            Upper bound of the parameter space
% Output: 
%   delta:         estimates of parameters (1 by 2*d)
%   fval:          the value of the objective function at the solution

% Compute the outcome for each configuration combination of x and delta

if strcmp(algorithm, 'LM-BFGS')
    % Construct the log likelihood function and score
    fcn = @(nuisance) ll(nuisance', x1, x2, data);
    grad = @(nuisance) score(nuisance', x1, x2, data);
    % Upper and lower bounds
    l  = lb;
    u  = ub;
    fun = @(nuisance) fminunc_wrapper(nuisance, fcn, grad);
    % Request very high accuracy for this test:
    opts2 = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
    opts2.printEvery = 100;
    opts2.x0 = delta_initial';
    % Run the algorithm:
    [delta, fval, ~] = lbfgsb(fun, l, u, opts2);
elseif strcmp(algorithm, 'BFGS')
    % Objective function and gradient
    f = @(delta_est) obj(delta_est, x1, x2, data);
    % Setting Options
    % Once gradient is supplied, switch to trust region algorithm
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
    %options2 = optimoptions('fminunc','Algorithm','quasi-newton');
    % Call fminunc
    [delta, fval] = fminunc(f, delta_initial, options);
end
end

function out = ll(delta, x1, x2, y)
% inputs:
% delta:  nuisance parameter; 1 by 6 (first 3 are for player 1)
% x1:     player 1 covariates (n by 3)
% x2:     player 2 covariates (n by 3)
% y:      entry-exit outcome for each market 
%
% output:
% out:    negative log likelihood

xdelta1 = x1*delta(1:3)'; % n by 1
xdelta2 = x2*delta(4:6)'; % n by 1
Phi1 = normcdf(xdelta1); 
Phi2 = normcdf(xdelta2);
% The following Q0 summarizes all possible combinations of entry-exit outcomes and covarites 
% The column corresponds to market outcomes 00, 01, 10 and 11
Q0 = [(1-Phi1).*(1-Phi2), (1-Phi1).*Phi2, Phi1.*(1-Phi2), Phi1.*Phi2];
% Now compute the negative log likelihood
mat_index = mkt_index(y);
out = -sum(log(Q0) .* mat_index, 'all');
end

function mat = mkt_index(y)
% create an n by 4 matrix that indexes which entry-exit outcome is realized
% for each market. Therefore each row sums up to 1. 

n = size(y, 1); % sample size
mat = zeros(n, 4);
mat(:, 1) = double((y==0));
mat(:, 2) = double((y==1));
mat(:, 3) = double((y==10));
mat(:, 4) = double((y==11));
end


% Compute z_delta and z_beta for each potential outcome
% Note: the code for the application is simplified because we do not
% compute power. Under the null, there is no subcases to consider. 
function z_delta = compute_z(delta, x1, x2)
% Inputs: 
%  delta:  nuisance parameter (1 by 6)
%          delta order: presence, size and constant
%
% Outputs: 
%  z_delta: n by 24 matrix with the following form
%                                LCC                            OA
%                      delta_1       delta_2   delta3   delta_4 ... delta_6
%           (0,0) (0,1) (1,0) (1,1)|         |        |
%   market1                        |         |        |
%   ...                            |         |        |
%   marketn                        |         |        |

n = size(x1, 1); % sample size
z_delta = zeros(n, 24); % Placeholders
beta = [0, 0];
% n by 2
xdelta = [x1*delta(1:3)', x2*delta(4:6)'];
phi = normpdf(xdelta);
Phi = normcdf(xdelta);
phi_beta = normpdf(xdelta + beta);
Phi_beta = normcdf(xdelta + beta);

% The following summarizes all possible combinations of entry-exit outcomes and covarites 
% The column corresponds to market outcomes 00, 01, 10 and 11
% Later will be dot mulitplied with an index matrix to compute the scores

% placeholders
zdelta1_10 = zeros(n, 3);
zdelta1_01 = zeros(n, 3);
zdelta2_10 = zeros(n, 3);
zdelta2_01 = zeros(n, 3);

% Event (0, 0) 
zdelta1_00 = -phi(:,1).*x1./(1-Phi(:,1)); % n by 3
zdelta2_00 = -phi(:,2).*x2./(1-Phi(:,2)); % n by 3
% Event (1, 1) 
zdelta1_11 = phi(:,1).*x1./Phi(:,1); % n by 3
zdelta2_11 = phi(:,2).*x2./Phi(:,2); % n by 3

var1 = Phi(:,1).*(1-Phi_beta(:,2));
z1 = Phi(:,1).*(1-Phi(:,2));
z2 = (1-Phi(:,1)).*Phi(:,2) + Phi(:,1) + Phi(:,2) - Phi(:,1).*Phi(:,2) - Phi_beta(:,1).*Phi_beta(:,2);
var2 = (z1.*z2 - Phi(:,2).*(1-Phi(:,1)).*z1)./(Phi(:,1) + Phi(:,2) - 2*Phi(:,1).*Phi(:,2));
var3 = Phi(:,1).*(1-Phi(:,2));
var4 = Phi_beta(:,1).*Phi(:,2);
var5 = Phi_beta(:,1).*Phi_beta(:,2);

for i = 1:n
    if var1(i) < var2(i) && var1(i) > var3(i) + var4(i) - var5(i) % (1, 0) chosen
        
        q1_01 = (1-Phi(i,1))*Phi(i,2) + Phi_beta(i,2)*(Phi(i,1)-Phi_beta(i,1));

        zdelta1_10(i,:) = x1(i,:).*phi(i,1)./Phi(i,1);

        zdelta2_10(i,:) = -x2(i,:).*phi_beta(i,2)/(1-Phi_beta(i,2));

        zdelta1_01(i,:) = (-x1(i,:)*Phi(i,2)*phi(i,1) + Phi_beta(i,2)*...
            x1(i,:)*(phi(i,1)-phi_beta(i,1)))./q1_01;

        zdelta2_01(i,:) = (x2(i,:)*(1-Phi(i,1))*phi(i,2) + x2(i,:)*...
            (Phi(i,1)-Phi_beta(i,1))*phi_beta(i,2))./q1_01;
        
    elseif var4(i) + var3(i) > var5(i) + var2(i) && var1(i) > var3(i) + var4(i) - var5(i)  % (0,1) chosen
        
        q1_10 = (1-Phi(i,2))*Phi(i,1) + Phi_beta(i,1)*(Phi(i,2)-Phi_beta(i,2));
        zdelta1_10(i,:) = x1(i,:).*( (1-Phi(i,2))*phi(i,1) + Phi(i,2)*phi_beta(i,1) - Phi_beta(i,2)*phi_beta(i,1) )...
            ./ q1_10;
        zdelta2_10(i,:) = x2(i,:) .* (-phi(i,2)*Phi_beta(i,1)+Phi(i,1)*phi(i,2)-...
            Phi_beta(i,1)*phi_beta(i,2))./q1_10;
        
        zdelta1_01(i,:) = -x1(i,:).*phi_beta(i,1)./(1-Phi_beta(i,1));
        
        zdelta2_01(i,:) = x2(i,:).*phi(i,2)./Phi(i,2);
        
    else % mixture
        denominator1 = Phi(i,1)+Phi(i,2)-Phi(i,1)*Phi(i,2)-Phi_beta(i,1)*Phi_beta(i,2);

        denominator2 = Phi(i,1)+Phi(i,2)-2*Phi(i,1)*Phi(i,2);
        
        zdelta1_10(i,:) = phi(i,1)*x1(i,:)/Phi(i,1) + (phi(i,1)*x1(i,:)...
            -phi(i,1)*x1(i,:)*Phi(i,2)-phi_beta(i,1)*Phi_beta(i,2)*...
            x1(i,:))/denominator1 - (phi(i,1)*x1(i,:)*(1-2*Phi(i,2)))/denominator2;
        
        zdelta2_10(i,:) = -phi(i,2)*x2(i,:)/(1-Phi(i,2)) + (phi(i,2)*...
            x2(i,:)-phi(i,2)*Phi(i,1)*x2(i,:)-phi_beta(i,2)*...
            Phi_beta(i,1)*x2(i,:))/denominator1 - (phi(i,2)*x2(i,:)*...
            (1-2*Phi(i,1)))/denominator2;


        zdelta1_01(i,:) = -phi(i,1)*x1(i,:)/(1-Phi(i,1)) + (phi(i,1)*...
            x1(i,:)-phi(i,1)*Phi(i,2)*x1(i,:)-phi_beta(i,1)*...
            Phi_beta(i,2)*x1(i,:))/denominator1 - (phi(i,1)*...
            x1(i,:)*(1-2*Phi(i,2)))/denominator2;

        zdelta2_01(i,:) = phi(i,2)*x2(i,:)/Phi(i,2) + (phi(i,2)*x2(i,:)-...
            phi(i,2)*Phi(i,1)*x2(i,:)-phi_beta(i,2)*Phi_beta(i,1)*...
            x2(i,:))/denominator1 - (phi(i,2)*x2(i,:)*(1-2*Phi(i,1)))/denominator2;
    end
end

z_delta(:, [1, 5, 9]) = zdelta1_00;
z_delta(:, [2, 6, 10]) = zdelta1_01;
z_delta(:, [3, 7, 11]) = zdelta1_10;
z_delta(:, [4, 8, 12]) = zdelta1_11;
z_delta(:, [13, 17, 21]) = zdelta2_00;
z_delta(:, [14, 18, 22]) = zdelta2_01;
z_delta(:, [15, 19, 23]) = zdelta2_10;
z_delta(:, [16, 20, 24]) = zdelta2_11;
end

% Score function
function g = score(delta, x1, x2, y)
zdelta = compute_z(delta, x1, x2); % n by 24
mat_index = mkt_index(y); % n by 4
ind = repmat(mat_index, [1, 6]) .* zdelta;
g = [sum(ind(:, 1:4), 'all'), sum(ind(:, 5:8), 'all'), sum(ind(:, 9:12), 'all'),...
    sum(ind(:, 13:16), 'all'), sum(ind(:, 17:20), 'all'), sum(ind(:, 21:24), 'all')];
% Since we add negative sign to use minimization algorithm, adjust the
% score accordingly
g = -g'; % 6 by 1
end

% for BFGS
function [f, g] = obj(delta, x1, x2, y)
f = ll(delta, x1, x2, y);
% gradient required
if nargout > 1 
    % Call compute_z_all to get zdelta, 2 by 16.
    g = score(delta, x1, x2, y);
end
end