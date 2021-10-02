% 10/1/2021 Shuowen Chen and Hiroaki Kaido
% Computes test statistic for application without discretization
function [Test, g_delta, varI, h, Vstar] = test_stat_no_discretization(delta_hat, data, x1, x2)
% inputs:
%   delta_hat:  first-step estimated delta
%   data:       market outcomes   
%   x1:         covariates of player 1, same definition for x2 (n by d)
% output:
%   Test:       test statstic
%   g_delta:    gn statistic
%   varI:       variance of gn statistic
%   h:          minimizer of quadratic program
%   Vstar:      minimum value of qudratic program
n = size(data, 1); % sample size

% 1. use compute_z (auxiliary fcn in this file) to get
% z_delta: n by 24 matrix that stores z for each delta under all 4 market outcome scenarios
% z_beta:  n by 8 matrix that stores z for each beta
[z_delta, z_beta] = compute_z(delta_hat', x1, x2);

% 2. use mlt_index (auxiliary fcn in this file) to get
% mat_index: n by 4 matrix, indicator of outcome realization for each
% market
mat_index = mkt_index(data);

% 3. compute the following individual score
%      p1      p2                p1                       p2
% [s_beta1, s_beta2]  [s_delta1 ... s_delta3]  [s_delta4 ... s_delta6]
s_beta = repmat(mat_index, [1, 2]) .* z_beta; 
s_delta = repmat(mat_index, [1, 6]) .* z_delta; 

% 4. compute the Cbeta and Cdelta
sum_s_beta = [sum(s_beta(:, 1:4), 'all'), sum(s_beta(:, 5:8), 'all')]; % 1 by 2
sum_s_delta = [sum(s_delta(:, 1:4), 'all'), sum(s_delta(:, 5:8), 'all'), ...
    sum(s_delta(:, 9:12), 'all'), sum(s_delta(:, 13:16), 'all'), ...
    sum(s_delta(:, 17:20), 'all'), sum(s_delta(:, 21:24), 'all')]; % 1 by 6

% Cdelta has three zero terms, leading to singularity issue of test
% statistic. This issue may be due to the fact that we use the RMLE whose 
% FOC tries to set C_delta to a vector close to 0.
% Therefore, we compute an analog of our test statistic while using only 
% the first term in g_n (i.e. using C_beta only)

each_score_beta = [sum(s_beta(:,1:4),2), sum(s_beta(:,5:8),2)];
varI = cov(each_score_beta); % covariance matrix
Cbeta = sum_s_beta/sqrt(n);
g_delta = Cbeta';
gtilde = varI^(-0.5)*g_delta;
Test = gtilde'*gtilde;
Vmin = @(x) quadform(varI, g_delta - x');
[h, Vstar] = fmincon(Vmin, [0,0], eye(2), zeros(2,1));
Test = Test - Vstar;
% 5. construct test statistic
% I_beta2 is 2 by 2; I_betadelta is 2 by 6, and I_delta2 is 6 by 6. 
%info_mat = [sum_s_beta'; sum_s_delta'] * [sum_s_beta, sum_s_delta];
%I_beta2 = info_mat(1:2, 1:2);
%I_delta2 = info_mat(3:end, 3:end);
%I_betadelta = info_mat(1:2, 3:end);

%I_delta2 = get_sigmatilde(I_delta2); % regularize by Andrews and Barwick (2012)

%g_delta = Cbeta' - I_betadelta/I_delta2 * Cdelta';
%varI = I_beta2 - I_betadelta/I_delta2 * I_betadelta';
%varI = get_sigmatilde(varI);

%gtilde = varI^(-0.5)*g_delta;
%Test = gtilde'*gtilde;
%Vmin = @(x) quadform(varI,g_delta-x);
%[~,Vstar] = fmincon(Vmin,[0;0],eye(2),zeros(2,1));
%Test = Test - Vstar;
end

% Auxiliary functions
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
function [z_delta, z_beta] = compute_z(delta, x1, x2)
% Inputs: 
%  delta:  nuisance parameter (2*3 by 1)
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
%
%  z_beta: a n by 8 matrix with the following form
%                       LCC               OA
%                      beta_1           belta_2  
%            (0,0) (0,1) (1,0) (1,1)|       
%   market1                         |
%   ...                             |
%   marketn                         |
% 
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
zbeta_10 = zeros(n, 2);
zbeta_01 = zeros(n, 2);

% Event (0, 0) 
zdelta1_00 = -phi(:,1).*x1./(1-Phi(:,1)); % n by 3
zdelta2_00 = -phi(:,2).*x2./(1-Phi(:,2)); % n by 3
zbeta_00 = zeros(n, 2);

% Event (1, 1) 
zdelta1_11 = phi(:,1).*x1./Phi(:,1); % n by 3
zdelta2_11 = phi(:,2).*x2./Phi(:,2); % n by 3
zbeta_11 = phi_beta./Phi_beta; % n by 2

var1 = Phi(:,1).*(1-Phi_beta(:,2));
z1 = Phi(:,1).*(1-Phi(:,2));
z2 = (1-Phi(:,1)).*Phi(:,2) + Phi(:,1) + Phi(:,2) - Phi(:,1).*Phi(:,2) - Phi_beta(:,1).*Phi_beta(:,2);
var2 = (z1.*z2 - Phi(:,2).*(1-Phi(:,1)).*z1)./(Phi(:,1) + Phi(:,2) - 2*Phi(:,1).*Phi(:,2));
var3 = Phi(:,1).*(1-Phi(:,2));
var4 = Phi_beta(:,1).*Phi(:,2);
var5 = Phi_beta(:,1).*Phi_beta(:,2);

for i = 1:n
    if var4(i) + var3(i) - var5(i) - var2(i) > 0 % (0,1) chosen
        q1_10 = (1-Phi(i,2))*Phi(i,1) + Phi_beta(i,1)*(Phi(i,2)-Phi_beta(i,2));
        zdelta1_10(i,:) = x1(i,:).*( (1-Phi(i,2))*phi(i,1) + Phi(i,2)*phi_beta(i,1) - Phi_beta(i,2)*phi_beta(i,1) )...
            ./ q1_10;
        zdelta2_10(i,:) = x2(i,:) .* (-phi(i,2)*Phi_beta(i,1)+Phi(i,1)*phi(i,2)-...
            Phi_beta(i,1)*phi_beta(i,2))./q1_10;
        
        zdelta1_01(i,:) = -x1(i,:).*phi_beta(i,1)./(1-Phi_beta(i,1));
        
        zdelta2_01(i,:) = x2(i,:).*phi(i,2)./Phi(i,2);
        
        zbeta_10(i, :) = [(Phi(i,2)-Phi_beta(i,2))*phi_beta(i,1)/q1_10,...
                          Phi_beta(i,1)*phi_beta(i,2)/q1_10];
        
        zbeta_01(i, :) = [-phi_beta(i, 1)/(1-Phi_beta(i, 1)), 0];
        
    elseif var1(i) - var2(i) < 0 % (1, 0) chosen
        q1_01 = (1-Phi(i,1))*Phi(i,2) + Phi_beta(i,2)*(Phi(i,1)-Phi_beta(i,1));

        zdelta1_10(i,:) = x1(i,:).*phi(i,1)./Phi(i,1);

        zdelta2_10(i,:) = -x2(i,:).*phi_beta(i,2)/(1-Phi_beta(i,2));

        zdelta1_01(i,:) = (-x1(i,:)*Phi(i,2)*phi(i,1) + Phi_beta(i,2)*...
            x1(i,:)*(phi(i,1)-phi_beta(i,1)))./q1_01;

        zdelta2_01(i,:) = (x2(i,:)*(1-Phi(i,1))*phi(i,2) + x2(i,:)*...
            (Phi(i,1)-Phi_beta(i,1))*phi_beta(i,2))./q1_01;
        
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
        
        zbeta_10(i, :) = [-phi_beta(i, 1)*Phi_beta(i, 2)/denominator1,...
                          -phi_beta(i, 2)*Phi_beta(i, 1)/denominator1];
        
        zbeta_01(i, :) = zbeta_10(i, :);
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

z_beta = zeros(n, 8);
z_beta(:, [1, 5]) = zbeta_00;
z_beta(:, [2, 6]) = zbeta_01;
z_beta(:, [3, 7]) = zbeta_10;
z_beta(:, [4, 8]) = zbeta_11;
end

% 3. Regularization using Andrews and Barwick
function SigmaTilde = get_sigmatilde(Sigma)
% Inputs
% Sigma: Covariance matrix
% Output
% SigmaTilde: Regularized covariance matrix based on Andrews and Barwick
sh_param = 0.05; 
sigma = sqrt(diag(Sigma));
Omega = Sigma./(sigma*sigma');
D = diag(diag(Sigma));
SigmaTilde = Sigma + max(sh_param-det(Omega), 0).*D;
end

% 4. Define quadratic form
function t = quadform(A, g)
gtilde = A^(-0.5)*g;
t = gtilde'*gtilde;
end