% 6/25/2020 Shuowen Chen and Hiroaki Kaido
% Computes the Tn test statistics in application
% The I_delta2 and varI are well behaved so that regularization is not
% needed
function Test = statapp(delta_hat, data, x1, x2)
% inputs:
%   delta_hat:  first-step estimated delta
%   data:       market outcomes   
%   x1:         covariates of player 1, same definition for x2 (n by d)
%   lambda:     regularization parameter 
% outputs:
%   Test: sup test statstics

[zdelta, zbeta] = compute_z_allapp([0, 0], delta_hat);
outcome = count(data, x1, x2);
Cbeta = zbeta * outcome / sqrt(size(data, 1));
Cdelta = zdelta * outcome / sqrt(size(data, 1));
% Compute the outer product approximation of information matrix and
% decomposed to four parts for orthogonalization
% I_beta2 is 2 by 2; I_betadelta and I_deltabeta are 2 by 6 (6 by 2),
% and I_delta2 is 6 by 6. 

% The following is intermediate for I_beta2. Dim is 4 by 64, where 64
% denotes all possible combinations of market outcome and covariates
betablock = [zbeta(1,:).*zbeta(1,:); zbeta(1,:).*zbeta(2,:);...
    zbeta(2,:).*zbeta(1,:); zbeta(2,:).*zbeta(2,:)];

% The following block is intermediate for I_delta2. Dim is 36 by 64.
% The first dimension is 36 since each player has 3 nuisance parameters and
% we need to consider each pairwise combination
deltablock = zeros(36, 64);
for i = 1:6
    for j = 1:6
        deltablock((j+6*(i-1)),:) = zdelta(i,:).*zdelta(j,:);
    end
end

% The following block is intermediate for I_deltabeta and I_betadelta (which
% is the transpose of the former). 
crossblock = zeros(12, 64);
for i = 1:2
    for j = 1:6
        crossblock((j+6*(i-1)),:) = zbeta(i,:).*zdelta(j,:);
    end
end

% Muliply by outcome of each combination
betablock = betablock * outcome / size(data, 1);
deltablock = deltablock * outcome / size(data, 1);
crossblock = crossblock * outcome / size(data, 1);

% Reshape
I_beta2 = reshape(betablock, 2, 2)'; % 2 by 2
I_delta2 = reshape(deltablock, 6, 6)'; % 6 by 6
I_betadelta = reshape(crossblock, 6, 2)'; % 2 by 6

% Check if I_delta2 is well behaved
% tol = 0.005;
% [~,D]=ldl(I_delta2);
% min(abs(diag(D))) > tol, hence no need to regularize

% Check if varI is well behaved
% [LI,DI] = ldl(varI);
% DI has no negative values, and min(abs(diag(DI))) > tol, no need for
% regularization

% Calculate the sup statistics
g_delta = Cbeta -  I_betadelta/I_delta2*Cdelta;
varI = I_beta2 - I_betadelta/I_delta2*I_betadelta';
Test = max(abs(varI^(-0.5)*g_delta));
end

% Auxiliary functions
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
% 2. Compute scores for each combination
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