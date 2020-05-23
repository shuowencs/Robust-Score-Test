% 5/21/2020 Shuowen Chen and Hiroaki Kaido
% Generates data 
% When beta is 0 (null), model is complete, no selection needed, use
% 'Complete' instead
% When beta is not 0 (alternative), model features multiple equilibria,
% use one of the three selection mechanisms to determine the outcome:
% (1) iid, (2) non iid (nonergodic), and (3) LFP.
function [data, x1, x2] = data_generation(beta, covariate, delta, n, S, selection)
%  Input:
%   beta:       structural parameters of interest (1 by 2)
%   covariate:  covariates (1 by 2) we draw from: (1, -1)
%   delta:      coefficients of covariates (nuisance parameters) (1 by 2)
%   n:          sample size
%   S:          number of simulations
%   selection:  equilibrium selection mechanism
%  Output:  
%   data:       the generated dataset (n by S). Each column denotes one
%               monte carlo simulation with sample size n. 
%   x1:         sampled x for the first player (n by S)
%   x2:         sampled x for the second player (n by S)
% set seed for replication
rng(123)
% generate a sequence of x 
x = randsample(covariate, n*2*S, true);
x1 = reshape(x(1:(n*S)), n, S);
x2 = reshape(x((n*S+1):length(x)), n, S);
% calculate xdelta
xdelta1 = x1*delta(1);
xdelta2 = x2*delta(2);

u1 = randn([n, S]);
u2 = randn([n, S]);

% Data generation under the null
if all(beta == 0) && strcmp(selection, 'Complete')
    data = NaN(n, S); % placeholder for data
    % Plug in the event values based on the model predictions
    data((u1 < -xdelta1) & (u2 < -xdelta2)) = 00;
    data((u1 > -xdelta1) & (u2 > -xdelta2)) = 11;
    data((u1 < -xdelta1) & (u2 > -xdelta2)) = 01;
    data((u1 > -xdelta1) & (u2 < -xdelta2)) = 10;
elseif all(beta == 0) && (strcmp(selection, 'Complete') == false)
    error('The model is complete under the null. Specify selection to be Complete.')
elseif ~all(beta == 0) && (strcmp(selection, 'Complete') == true)
    error('The model is incomplete given parameter specifications. Use IID, Non IID or LFP in selection')
end
% Data generation under the alternative
if strcmp(selection, 'IID')
    % iid selection mechanism, selects (1,0) out of {(1,0),(0,1)} if an
    % iid Bernoulli r.v vi takes 1, and (0,1) if vi takes 0.
    % Random numbers from standard normal distribution
    % Define the Bernoulli r.v: threshold value to create Bernoulli variables
    alpha = 0.05; 
    vi = rand(n,S) <= alpha;
    % placeholder for data
    G_result = NaN(n, S);
    % Plug in the event values based on the model predictions
    G_result((u1 < -xdelta1) & (u2 < -xdelta2)) = 00;
    G_result((u1 > -xdelta1-beta(1)) & (u2 > -xdelta2-beta(2))) = 11;
    G_result((u1 < -xdelta1) & (u2 > -xdelta2)) = 01;
    G_result((-xdelta1 < u1) & (u1 < -xdelta1-beta(1)) & ...
        (u2 > -xdelta2-beta(2))) = 01;
    G_result((u1 > -xdelta1-beta(1)) & (u2 < -xdelta2-beta(2))) = 10;
    G_result((-xdelta1 < u1) & (u1 < -xdelta1-beta(1)) & ...
        (u2 < -xdelta2)) = 10;
    % Use vi to select outcome from multiple equilibria
    G_result((-xdelta1 < u1) & (u1 < -xdelta1-beta(1)) & ...
        (-xdelta2 < u2) & (u2 < -xdelta2-beta(2)) & vi) = 10;
    G_result((-xdelta1 < u1) & (u1 < -xdelta1-beta(1)) & ...
        (-xdelta2 < u2) & (u2 < xdelta2-beta(2)) & (~vi)) = 01;
    data = G_result;
end

if strcmp(selection,'Non IID')
    % Define N^{*}_{k} as an increasing sequence of integers
    get_hi = @(i) 2.^(2.^ceil(log(log(i)./log(2))./log(2)));
    nL = get_hi(n);
    % Total Data Set for u, x and xdelta
    % NOTE: in Non IID case, we actually work with data larger than the
    % specified sample size due to the fact that the last cluster might
    % contain a number of obs such that the total obs of all cluster
    % exceeds the sample size n. At the end of the DGP process we subset
    % sample size n. The economic interpretation is that we only observe a
    % faction of the obs in the last cluster. As a result, we define larger
    % u, x, and hence xdelta.
    u1L = randn([nL,S]);
    u2L = randn([nL,S]);
    xL = randsample(covariate, nL*2*S, true);
    x1L = reshape(xL(1:(nL*S)), nL, S);
    x2L = reshape(xL((nL*S+1):length(xL)), nL, S);
    xLdelta1 = x1L*delta(1);
    xLdelta2 = x2L*delta(2);
    hi_previous = 1;
    datL = zeros(nL, S);
    x_col = [reshape(x1L, nL*S, 1), reshape(x2L, nL*S, 1)];
    % Reshape the x_col matrix to be of dim nL by 2*S, where the first
    % two columns are x1 and x2 for first MC, etc. This will be used for
    % counting the # of each of the four configs of X.
    x_nS = NaN(nL, 2*S);
    for i = 1:S
        x_nS(:, (2*i-1):(2*i)) = x_col((1+(i-1)*nL):(i*nL), :);
    end
    % The loop won't break until the total number of obs exceeds n
    while 1
        % the last obs in the current cluster
        hi_current = get_hi(hi_previous+1);
        if hi_current == hi_previous
            error('K error')
        end
        % Count non-incomplete area number of (1,0) and (0,1), consider
        % unique prediction area
        % The number of observations in the current cluster
        cluster_current = (hi_previous+1):hi_current;
        
        % Count # of (1,0)
        % NOTE: we start counting from the first until the current cluster
        G10_1 = sum(u1L(1:hi_current,:)>-xLdelta1(1:hi_current,:) & ...
            u1L(1:hi_current,:)<-xLdelta1(1:hi_current,:)-beta(1) & ...
            u2L(1:hi_current,:)<-xLdelta2(1:hi_current,:));
        G10_2 = sum(u1L(1:hi_current,:)>-xLdelta1(1:hi_current,:)-beta(1) & ...
            u2L(1:hi_current,:)<-xLdelta2(1:hi_current,:)-beta(2));
        % Two parts in correspondence give us (1,0) outcome
        G10cum = G10_1 + G10_2;
        
        % Count # of (0,1)
        G01_1 = sum(u1L(1:hi_current,:)>-xLdelta1(1:hi_current,:) & ...
            u1L(1:hi_current,:)<-xLdelta1(1:hi_current,:)-beta(1) & ...
            u2L(1:hi_current,:)>-xLdelta2(1:hi_current,:)-beta(2));
        G01_2 = sum(u1L(1:hi_current,:)<-xLdelta1(1:hi_current,:) & ...
            u2L(1:hi_current,:)>-xLdelta2(1:hi_current,:));
        % Two parts in correspondence give us (0,1) outcome
        G01cum = G01_1 + G01_2;
        
        % Calculate the ratio of frequencies
        ratio = G10cum./(G10cum + G01cum);
        
        % Calcualte the threshold Lambda
        % Now we count the # of each of the four possible configs of X in
        % each MC from the 1st obs until the last obs of the current cluster
        x_current = x_nS(1:hi_current, :);
        % placeholder
        config_count = NaN(4, S);
        for i = 1:S
            config_count(1, i) = sum(x_current(:, 2*i-1)==1 & x_current(:, 2*i)==1);
            config_count(2, i) = sum(x_current(:, 2*i-1)==1 & x_current(:, 2*i)==-1);
            config_count(3, i) = sum(x_current(:, 2*i-1)==-1 & x_current(:, 2*i)==1);
            config_count(4, i) = sum(x_current(:, 2*i-1)==-1 & x_current(:, 2*i)==-1);
        end
        
        % Now calculate the lower probabilities of events (1,0) and (0,1),
        % conditional on configuration of X. This contains 8 analytical
        % expressions. Call function lowerprob.
        X = [1, 1; 1, -1; -1, 1; -1, -1];
        % nu_10 and nu_01 are both 4 by 1
        [nu_10, nu_01] = lowerprob(X, delta, beta);
        % Now we calculate the threshold Lambda_hi (1 by S matrix)
        % The numerator
        numerator_hi = sum(config_count.*nu_10, 1);
        % The denominator
        denominator_hi = numerator_hi + sum(config_count.*nu_01, 1);
        Lambda_hi = numerator_hi./denominator_hi;
        % Calculate the ratio (1,0) over (1,0) and (0,1)
        vii= ratio > Lambda_hi;
        % Need to repeat the dummy for each observation within the
        % cluster
        vi=repmat(vii,length(cluster_current),1);
        
        % Create logic variable
        G00 = (u1L(cluster_current,:)<-xLdelta1(cluster_current,:) & ...
            u2L(cluster_current,:)<-xLdelta2(cluster_current,:));
        
        G11 = (u1L(cluster_current,:)>-xLdelta1(cluster_current,:)-beta(1) & ...
            u2L(cluster_current,:)>-xLdelta2(cluster_current,:)-beta(2));
        
        G10_1 = (u1L(cluster_current,:)>-xLdelta1(cluster_current,:) & ...
            u1L(cluster_current,:)<-xLdelta1(cluster_current,:)-beta(1) &...
            u2L(cluster_current,:)<-xLdelta2(cluster_current,:));
        
        G10_2 = (u1L(cluster_current,:)>-xLdelta1(cluster_current,:)-beta(1) &...
            u2L(cluster_current,:)<-xLdelta2(cluster_current,:)-beta(2));
        
        G01_1 = (u1L(cluster_current,:)>-xLdelta1(cluster_current,:) & ...
            u1L(cluster_current,:)<-xLdelta1(cluster_current,:)-beta(1)...
            & u2L(cluster_current,:)>-xLdelta2(cluster_current,:)-beta(2));
        
        G01_2 = (u1L(cluster_current,:)<-xLdelta1(cluster_current,:) & ...
            u2L(cluster_current,:)>-xLdelta2(cluster_current,:));
        
        % Generate cluster data
        dat_cluster = NaN(size(G11));
        dat_cluster(G11) = 11;
        dat_cluster(G00) = 00;
        dat_cluster(G10_1) = 10;
        dat_cluster(G10_2) = 10;
        dat_cluster(G01_1) = 01;
        dat_cluster(G01_2) = 01;
        
        % select outcome from Bernoulli vi
        dat_cluster((~G11) & (~G00) & (~G10_1) & (~G10_2) & (~G01_1) & (~G01_2) & vi) = 10;
        dat_cluster((~G11) & (~G00) & (~G10_1) & (~G10_2) & (~G01_1) & (~G01_2) & (~vi)) = 01;
        datL(cluster_current,:) = dat_cluster;
        
        % Update the hi
        hi_previous = hi_current;
        
        % If the total number of observations is bigger than n, break
        % the while loop
        if hi_current >= n; break; end
    end
    % subset the first n obs for each MC
    data = datL(1:n,:);
    x1 = x1L(1:n, :);
    x2 = x2L(1:n, :);
end

if strcmp(selection, 'LFP')
    % generate data from a least favorable distribution
    % Draws an outcome sequence from Q_0 if beta = [0,0] and Q_1 o.w.
    data = NaN(n*S, 1);
    if all(beta == 0)
        % gen Q0 distributions for each obs
        [Q, ~] = get_LFP(beta, xdelta1, xdelta2);
    elseif all(beta < 0)
        % gen Q1 distributions for each obs
        [~, Q] = get_LFP(beta, xdelta1, xdelta2);
    end
    % cumulative probability from the Q dist for each observation.
    cdf_Q = cumsum(Q, 2);
    % Generate uniform distribution (probabilities) matrix
    unif_prob = rand(n*S, 1);
    % Generate iid samples from Q1
    data((unif_prob >= 0) & (unif_prob < cdf_Q(:,1))) = 00;
    data((unif_prob >= cdf_Q(:,1)) & (unif_prob < cdf_Q(:,2))) = 11;
    data((unif_prob >= cdf_Q(:,2)) & (unif_prob < cdf_Q(:,3))) = 10;
    data((unif_prob >= cdf_Q(:,3)) & (unif_prob < cdf_Q(:,4))) = 01;
    data = reshape(data, n, S);
end
end

% Define two auxiliary functions.
function [Q0, Q1] = get_LFP(beta, xdelta1, xdelta2)
% Purpose: The function calculates the analytical LFP given beta and Xdelta. 
% Inputs:
%   beta:  structural parameters of interest
%   xdelta1: n by S matrix
%   xdelta2: n by S matrix
% Outputs:
%   Q0:     Distribution under the null [q0_00, q0_11, q0_10, q0_01],
%           Dimension: nS by 4. 
%   Q1:     Distribution under the alternative [q1_00, q1_11, q1_10, q1_01],
%           Dimension: nS by 4. 

dim = size(xdelta1);
% sample size
n = dim(1);
% number of MC
S = dim(2);
% Each of the following four objects are n by S. 
Phi1 = normcdf(xdelta1);
Phi2 = normcdf(xdelta2);
Phi_beta1 = normcdf(xdelta1 + beta(1));
Phi_beta2 = normcdf(xdelta2 + beta(2));
% Within a MC, Q0 is identical arocss all subcases. 
% The Q0 below is a n*S by 4 matrix 
Q0 = [reshape((1-Phi1).*(1-Phi2),n*S,1), reshape(Phi1.*Phi2,n*S,1),...
    reshape(Phi1.*(1-Phi2),n*S,1), reshape((1-Phi1).*Phi2,n*S,1)];
% Depending on the values of beta, we have three cases to
% consider. We have six conditions (two of which are repetitive) to
% consider the regions. Given each local alternative, we can evaluate
% the conditions and thus determine which subregion the given local
% alternative belongs to. First define the 5 variables used in the
% conditions. All 5 variables are n by S. 
var1 = Phi1.*(1-Phi_beta2);
z1 = Phi1.*(1-Phi2);
z2 = Phi2.*(1-Phi1) + Phi1 + Phi2 - Phi1.*Phi2 - Phi_beta1.*Phi_beta2;
var2 = reshape((z1.*z2 - Phi2.*(1-Phi1).*z1)./(Phi1 + Phi2 - 2*Phi1.*Phi2), n*S, 1);
var3 = reshape(Phi1.*(1-Phi2), n*S, 1);
var4 = reshape(Phi_beta1.*Phi2, n*S, 1);
var5 = reshape(Phi_beta1.*Phi_beta2, n*S, 1);
var1 = reshape(var1, n*S, 1);
% for condition 1, instead of using ==, we need to take into consideration
% the floating-point arithmetic using a tolerance
tol = eps(1);
% Case1: Mixture over (1,0) and (0,1)
cond1 = (var1 - var2 > 0 | abs(var1 - var2) <= tol) & ...
    (var2 - var3 - var4 + var5 > 0 | abs(var2 - var3 - var4 + var5) <= tol);
% Case2: (0,1) chosen
cond2 = var4 + var3 - var5 - var2 > 0 & var1 - var3 - var4 + var5 > 0;
% Case 3: (1,0) chosen
cond3 = var1 - var2< 0 & var1 - var3 - var4 + var5 > 0;
% stack up the evaluations
evaluation = [cond1, cond2, cond3];
% Now based on the evaluation, calculate Q1 for each obs. Instead of
% evaluating each obs case-by-case (using if else), we first calculate Q1
% in each of the three subcases. This gives us three nS by 4 matrices. Then
% we dot multiply each of them with the corresponding column in the evaluation 
% matrix and sum them up. This gives us the nS by 4 matrix of LFP for all obs 
% (S Monte Carlos each with sample size n).

% First calculate Q1 in each case
% Q1 in case 1: Mixture over (1, 0) and (0, 1)
Q1_1 = [reshape((1-Phi1).*(1-Phi2), n*S, 1),...
        reshape(Phi_beta1.*Phi_beta2, n*S, 1),...
        reshape((z1.*z2-Phi2.*(1-Phi1).*z1)./(Phi1+Phi2-2*Phi1.*Phi2), n*S, 1),...
        reshape(Phi1+Phi2-Phi1.*Phi2-Phi_beta1.*Phi_beta2 - ...
        (z1.*z2-Phi2.*(1-Phi1).*z1)./(Phi1+Phi2-2*Phi1.*Phi2), n*S, 1)];
% Q1 in case 2: (0, 1) chosen
Q1_2 = [reshape((1-Phi1).*(1-Phi2), n*S, 1),...
        reshape(Phi_beta1.*Phi_beta2, n*S, 1), ...
        reshape((1-Phi2).*Phi1+Phi_beta1.*(Phi2-Phi_beta2), n*S, 1),...
        reshape((1-Phi_beta1).*Phi2, n*S, 1)];
% Q1 in case 3: (1, 0) chosen
Q1_3 = [reshape((1-Phi1).*(1-Phi2), n*S, 1),...
      reshape(Phi_beta1.*Phi_beta2, n*S, 1),...
      reshape(Phi1.*(1-Phi_beta2), n*S, 1),...
      reshape((1-Phi1).*Phi2+Phi_beta2.*(Phi1-Phi_beta1), n*S, 1)];
% Now get Q1 for each obs (nS by 4)
Q1 = evaluation(:,1).*Q1_1 + evaluation(:,2).*Q1_2 + evaluation(:,3).*Q1_3;
end

function [nu_10, nu_01] = lowerprob(X, delta, beta)
% Purpose: calculate the (analytical) lower probabilities of events (1, 0)
% and (0, 1) conditional on X. Used in non IID selection 
% 
% Inputs: 
%   X: covariates [x1, x2]. Should be [1, 1; 1, -1; -1, 1; -1, -1]
%   delta: coefficients of X [delta1, delta2]
%   beta: structural parameters of interest [beta1, beta2]
% Outputs: 
%   nu_10: lower probabilities of event (1, 0)
%   nu_01: lower probabilities of event (0, 1)
%   NOTE: both outputs are 4 by 1

Xdelta = X.*delta;
nu_10 = (1-normcdf(Xdelta(:, 2))).*normcdf(Xdelta(:, 1)) + ...
    normcdf(Xdelta(:, 1)+beta(1)).*(normcdf(Xdelta(:, 2)) - normcdf(Xdelta(:, 2)+beta(2)));
nu_01 = (1-normcdf(Xdelta(:, 1))).*normcdf(Xdelta(:, 2)) + ...
    normcdf(Xdelta(:, 2)+beta(2)).*(normcdf(Xdelta(:, 1)) - normcdf(Xdelta(:, 1)+beta(1)));
end