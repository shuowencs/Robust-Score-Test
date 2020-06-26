% 6/16/2020 Shuowen Chen and Hiroaki Kaido
% Testing function for application (using L-BFGS-M and BFGS)
function [delta, fval] = rmleapp(delta_initial, algorithm, data, x1, x2, lb, ub)
% Input:
%   delta_initial: Initial guess of delta
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
outcome = count(data, x1, x2);
if strcmp(algorithm, 'LM-BFGS')
    % Construct the log likelihood function and score
    fcn = @(nuisance) loglikelihood(nuisance, outcome);
    grad = @(nuisance) score(nuisance, outcome);
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
    f = @(delta_est) obj(delta_est, outcome);
    % Setting Options
    % Once gradient is supplied, switch to trust region algorithm
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
    % Call fminunc
    [delta, fval] = fminunc(f, delta_initial', options);
end
end

% ---- Auxiliary functions ---- % 

% 1. Count the frequencies of combination of outcome and covariates
function output = count(data, x1, x2)
% turn into binary code 
% note we exclude the constants here since they always show up
temp = 10000*data + 1000*x1(:,1) + 100*x1(:,2) + 10*x2(:,1) + x2(:,2);
% turn into decimal number (0-63) so as to do counting using for loop
temp = bin2dec(num2str(temp));
% placeholder
output = NaN(64, 1);
for i = 0 : 63
    output(i+1) = sum(temp == i);
end
end

% 2. Compute Least favorable pair under the null
function Q0 = get_LFP(delta)
% Purpose: The function calculates the analytical LFP given and Xdelta. 
% Inputs:
%   x1:     covariate of player 1; n by 3 matrix
%   x2:     covariate of player 2; n by 3 matrix
%   delta:  nuisance parameter; 6 by 1 (first 3 are for player 1)
% Outputs:
%   Q0:     Distribution under the null [q0_00, q0_01, q0_10, q0_11],
%           Dimension: 16 by 4. Each market outcome admits four possibilities 
%           as there are four combinations of covariates
x = [0,0,1; 0,1,1; 1,0,1; 1,1,1];
xdelta1 = x*delta(1:3);
xdelta2 = x*delta(4:6);
% Each of the following four objects are 16 by 1. 
Phi1 = repelem(normcdf(xdelta1), 4, 1);
Phi2 = repmat(normcdf(xdelta2), [4, 1]);
% The Q0 below is a 16 by 4 matrix 
Q0 = [(1-Phi1).*(1-Phi2), (1-Phi1).*Phi2, Phi1.*(1-Phi2), Phi1.*Phi2];
end

% 3. Negative Log likelihood function
function f = loglikelihood(delta, outcome)
% delta should be 6 by 1
q0 = reshape(get_LFP(delta), 64, 1);
% objective function
f = -sum(outcome .* log(q0));
end

% 4. scores for each single market
function [z_delta, z_beta] = compute_z_allapp(beta, delta)
% Inputs: 
%  beta:   structural parameter of interest (1 by 2 vector)
%  delta:  nuisance parameter (2*3 by 1)
%          delta order: presence, size and constant
%
% Outputs: 
%  z_delta: a 6 by 64 matrix with the following form
%   outcome                  (0,0)                            (0,1) (1,0) (1,1)
%   covariate (0,0,1)(0,0,1) (0,0,1)(0,1,1) ...(1,1,1)(1,1,1)|     |     |     |
%   delta1_1                                                 |     |     |     |
%   delta1_2                                                 |     |     |     |
%   delta1_3                                                 |     |     |     |
%   delta2_1                                                 |     |     |     |
%   delta2_2                                                 |     |     |     |
%   delta2_3                                                 |     |     |     |
%
%  z_beta: a 2 by 64 matrix with the following form
%   outcome                  (0,0)                            (0,1) (1,0) (1,1)
%   covariate (0,0,1)(0,0,1) (0,0,1)(0,1,1) ...(1,1,1)(1,1,1)|     |     |     |
%   beta1_1                                                  |     |     |     |
%   beta1_2                                                  |     |     |     |
% 
% Placeholders:
z_delta = zeros(6, 64);
z_beta = zeros(2, 64);
% Construct the 4 possible covariate configurations
% Note: each player has the same 4 possible configurations, but when
% interacted there will be a totality of 16 cases
x = [0,0,1;0,1,1;1,0,1;1,1,1]; % constant is the 3rd column
% 4 possible combinations for x*delta for two players (4 by 2)
xdelta = [x*delta(1:3), x*delta(4:6)];
phi = normpdf(xdelta);
Phi = normcdf(xdelta);
phi_beta = normpdf(xdelta + beta);
Phi_beta = normcdf(xdelta + beta);

% For events (0, 0) and (1, 1) no subcase consideration
% zbeta_11 is 4 by 2
zbeta_11 = phi_beta./Phi_beta;
zbeta_00 = zeros(2, 16);
% Scores related to delta now features 4 possible combinations of 
% covariates and 3 coefficients (each is 4 by 3)
zdelta1_00 = -phi(:,1).*x./(1-Phi(:,1));
zdelta2_00 = -phi(:,2).*x./(1-Phi(:,2));
zdelta1_11 = phi_beta(:,1).*x./Phi_beta(:,1);
zdelta2_11 = phi_beta(:,2).*x./Phi_beta(:,2);
% Need to change Phi, phi and x to accommodate the combined considerations
Phi1 = repelem(Phi(:,1),4);
Phi2 = repmat(Phi(:,2),[4,1]);
phi1 = repelem(phi(:,1),4);
phi2 = repmat(phi(:,2),[4,1]);
Phi_beta1 = repelem(Phi_beta(:,1),4);
Phi_beta2 = repmat(Phi_beta(:,2),[4,1]);
phi_beta1 = repelem(phi_beta(:,1),4);
phi_beta2 = repmat(phi_beta(:,2),[4,1]);
x_first = repelem(x,4,1);
x_second = repmat(x,[4,1]);
% For events (1, 0) and (0, 1), for each configuration of Xdelta, 
% there are three subcases to consider. We have six conditions 
% (two of which are repetitive) to consider the regions. 
% Given each local alternative, we can evaluate the conditions and thus
% determine which subregion the given local alternative belongs to.
% First define the 5 variables used in the conditions
% Each is 16 by 1
var1 = Phi1.*(1-Phi_beta2);
z1 = Phi1.*(1-Phi2);
z2 = (1-Phi1).*Phi2 + Phi1 + Phi2 - Phi1.*Phi2 - Phi_beta1.*Phi_beta2;
var2 = (z1.*z2 - Phi2.*(1-Phi1).*z1)./(Phi1 + Phi2 - 2*Phi1.*Phi2);
var3 = Phi1.*(1-Phi2);
var4 = Phi_beta1.*Phi2;
var5 = Phi_beta1.*Phi_beta2;

% Pre-allocation. Note: each player's zdelta has 3 components
zdelta_0110 = zeros(6, 32);
zbeta_0110 = zeros(2, 32);

for i = 1:16 % loop over combined covariate configurations
    % (0,1) chosen
    if var4(i) + var3(i) - var5(i) > var2(i) && var1(i) > var3(i) + var4(i) - var5(i) 

        q1_10 = Phi1(i)*(1-Phi2(i))+Phi_beta1(i)*Phi2(i)-Phi_beta1(i)*Phi_beta2(i);

        zbeta1_10 = (Phi2(i)*phi_beta1(i)-Phi_beta2(i)*phi_beta1(i))/q1_10;

        zbeta2_10 = -Phi_beta1(i)*phi_beta2(i)/q1_10;

        zdelta1_10 = x_first(i,:) * (phi1(i)*(1-Phi2(i))+Phi2(i)*phi_beta1(i)-...
            Phi_beta2(i)*phi_beta1(i))/q1_10;

        zdelta2_10 = x_second(i,:) * (-phi2(i)*Phi_beta1(i)+Phi1(i)*phi2(i)-...
            Phi_beta1(i)*phi_beta2(i))/q1_10;

        zdelta1_01 = -x_first(i,:)*phi_beta1(i)/(1-Phi_beta1(i));
        
        zdelta2_01 = x_second(i,:)*phi2(i)/Phi2(i);

        zbeta1_01 = -phi_beta1(i)/(1-Phi_beta1(i));

        zbeta2_01 = 0; 
    % (1,0) chosen    
    elseif var1(i) < var2(i) && var1(i) - var3(i) - var4(i) + var5(i) > 0 

        q1_01 = (1-Phi1(i))*Phi2(i) + Phi_beta2(i)*(Phi1(i)-Phi_beta1(i));

        zdelta1_10 = x_first(i,:)*phi1(i)/Phi1(i);

        zdelta2_10 = -x_second(i,:)*phi_beta2(i)/(1-Phi_beta2(i));

        zbeta1_10 = 0;

        zbeta2_10 = -phi_beta2(i)/(1-Phi_beta2(i));

        zdelta1_01 = (-x_first(i,:)*Phi2(i)*phi1(i) + Phi_beta2(i)*...
            x_first(i,:)*(phi1(i)-phi_beta1(i)))/q1_01;

        zdelta2_01 = (x_second(i,:)*(1-Phi1(i))*phi2(i) + x_second(i,:)*...
            (Phi1(i)-Phi_beta1(i))*phi_beta2(i))/q1_01;

        zbeta1_01 = -Phi_beta2(i)*phi_beta1(i)/q1_01;

        zbeta2_01 = (Phi1(i)-Phi_beta1(i))*phi_beta2(i)/q1_01;
    else
        % Mixture over (1,0) and (0,1)
        denominator1 = Phi1(i)+Phi2(i)-Phi1(i)*Phi2(i)-Phi_beta1(i)*Phi_beta2(i);

        denominator2 = Phi1(i)+Phi2(i)-2*Phi1(i)*Phi2(i);
        
        zdelta1_10 = phi1(i)*x_first(i,:)/Phi1(i) + (phi1(i)*x_first(i,:)...
            -phi1(i)*x_first(i,:)*Phi2(i)-phi_beta1(i)*Phi_beta2(i)*...
            x_first(i,:))/denominator1 - (phi1(i)*x_first(i,:)*(1-2*Phi2(i)))/denominator2;
        
        zdelta2_10 = -phi2(i)*x_second(i,:)/(1-Phi2(i)) + (phi2(i)*...
            x_second(i,:)-phi2(i)*Phi1(i)*x_second(i,:)-phi_beta2(i)*...
            Phi_beta1(i)*x_second(i,:))/denominator1 - (phi2(i)*x_second(i,:)*...
            (1-2*Phi1(i)))/denominator2;

        zbeta1_10 = -phi_beta1(i)*Phi_beta2(i)/denominator1;

        zbeta2_10 = -Phi_beta1(i)*phi_beta2(i)/denominator1;

        zdelta1_01 = -phi1(i)*x_first(i,:)/(1-Phi1(i)) + (phi1(i)*...
            x_first(i,:)-phi1(i)*Phi2(i)*x_first(i,:)-phi_beta1(i)*...
            Phi_beta2(i)*x_first(i,:))/denominator1 - (phi1(i)*...
            x_first(i,:)*(1-2*Phi2(i)))/denominator2;

        zdelta2_01 = phi2(i)*x_second(i,:)/Phi2(i) + (phi2(i)*x_second(i,:)-...
            phi2(i)*Phi1(i)*x_second(i,:)-phi_beta2(i)*Phi_beta1(i)*...
            x_second(i,:))/denominator1 - (phi2(i)*x_second(i,:)*(1-2*Phi1(i)))/denominator2;

        zbeta1_01 = -phi_beta1(i)*Phi_beta2(i)/denominator1;

        zbeta2_01 = -phi_beta2(i)*Phi_beta1(i)/denominator1;
    end
    % stacking up the results
    zbeta_0110(:, 16+i) = [zbeta1_01; zbeta2_01];
    zbeta_0110(:, 32+i) = [zbeta1_10; zbeta2_10];
    zdelta_0110(:, 16+i) = [zdelta1_01'; zdelta2_01'];
    zdelta_0110(:, 32+i) = [zdelta1_10'; zdelta2_10'];
end

% Fill in the rest of z_beta
z_beta(:,1:16) = zbeta_00;
z_beta(:,49:end) = repelem(zbeta_11', 1, 4);
% Fill in the rest of z_delta
z_delta(1:3,1:16) = repelem(zdelta1_00', 1, 4);
z_delta(4:end,1:16) = repelem(zdelta2_00', 1, 4);
z_delta(1:3,49:end) = repelem(zdelta1_11', 1, 4);
z_delta(4:end,49:end) = repelem(zdelta2_11', 1, 4);
end

% 5. Score function
function g = score(delta, outcome)
% delta should be 2 by 1
% Call compute_z_all to get zdelta, 2 by 16.
[zdelta, ~] = compute_z_allapp([0, 0], delta);
% Compute sum of score functions (2 by 1)
g = zdelta * outcome;
% Since we add negative sign to use minimization algorithm, adjust the
% score accordingly
g = -g;
end

% 6 For BFGS: objective function and corresponding gradient
function [f, g] = obj(delta, outcome)
% compute the PMF for each market given covariates and delta (n by 4):
% [q0_00, q0_11, q0_10, q0_01]
% delta should be 6 by 1
q0 = reshape(get_LFP(delta), 64, 1);
% objective function
f = -sum(outcome .* log(q0));
% gradient required
if nargout > 1 
    % Call compute_z_all to get zdelta, 2 by 16.
    [zdelta, ~] = compute_z_allapp([0, 0], delta);
    % Compute sum of score functions (2 by 1)
    g = zdelta * outcome;
    % Since we add negative sign to use minimization algorithm, adjust the
    % score accordingly
    g = -g;
end
end