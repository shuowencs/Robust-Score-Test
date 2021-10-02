% 5/7/2021 Shuowen Chen and Hiroaki Kaido
% Acknowledgement: based on code by Bugni, Canay and Shi (2017, QE)
% Data manipulation to produce sample moment conditions for the game in
% Chen and Kaido (2021)

% We consider two types of moments: 
% "agg":    aggregates the covariates across markets 4 moment (in)equalities
% "disagg": conditional on the value of covariates, has 16 ( because the 
% support of x is {(-1,-1), (-1,1), (1,-1), (1,1)} )

function [mbar_std, mData, xi] = dataOp(delta, beta, data, kappa, X1, X2, mtype)
%{
 -Inputs-
 delta:    nuisance parameter
 beta:     the parameter of interest
 data:     dataset (n by S)
 kappa:    tuning parameter
 X1:       covariate for player 1 (n by S)
 X2:       covariate for player 2 (n by S)
 mtype:    whether we use aggregate or disaggregate moments

 -Outputs- 
 Let k denote the total number of moments
 mbar_std: standardized sample moment inequalities and equalities (k by S)
 mData:    data for sample moments (n by k by S). It is used for
           constructing objective function for DR and PR tests.  
 xi:       slackness measure in GMS (k by S)
%}

% sample size (number of markets);
n = size(data, 1);
% number of simulations
S = 1;
% S = size(data, 2);

% Model predictions (n by S)
% (1) equalities
modelP00 = normcdf(-X1*delta(1)) .* normcdf(-X2*delta(2));
modelP11 = normcdf(X1*delta(1)+beta(1)) .* normcdf(X2*delta(2)+beta(2));
% (2) inequalities
% region of multiplicity
mul = (normcdf(-X1*delta(1)-beta(1)) - normcdf(-X1*delta(1))) .* (normcdf(-X2*delta(2)-beta(2)) - normcdf(-X2*delta(2)));
modelP10_ub = normcdf(X1*delta(1)) .* normcdf(-X2*delta(2)-beta(2));
modelP10_lb = modelP10_ub - mul;

if strcmp(mtype, 'disagg')
    % Construct aggregate sample moments
    dataP00x1 = zeros(n, S); % placeholders for outcome and covariate combo
    dataP00x2 = zeros(n, S);
    dataP00x3 = zeros(n, S);
    dataP00x4 = zeros(n, S); 
    dataP11x1 = zeros(n, S);
    dataP11x2 = zeros(n, S);
    dataP11x3 = zeros(n, S);
    dataP11x4 = zeros(n, S);
    % outcome (1,0) will be used twice for the inequalities
    dataP10x1 = zeros(n, S); 
    dataP10x2 = zeros(n, S);
    dataP10x3 = zeros(n, S);
    dataP10x4 = zeros(n, S);
    for s = 1:S % loop over each simulation
        % From the data, count frequencies of possible outcomes (0, 0), (1, 1)
        % and (1, 0)
        dataP00x1(data(:, s) == 0 & X1(:, s) == 1 & X2(:, s) == 1, s) = 1;
        dataP00x2(data(:, s) == 0 & X1(:, s) == 1 & X2(:, s) == -1, s) = 1;
        dataP00x3(data(:, s) == 0 & X1(:, s) == -1 & X2(:, s) == 1, s) = 1;
        dataP00x4(data(:, s) == 0 & X1(:, s) == -1 & X2(:, s) == -1, s) = 1;
        dataP11x1(data(:, s) == 11 & X1(:, s) == 1 & X2(:, s) == 1, s) = 1;
        dataP11x2(data(:, s) == 11 & X1(:, s) == 1 & X2(:, s) == -1, s) = 1;
        dataP11x3(data(:, s) == 11 & X1(:, s) == -1 & X2(:, s) == 1, s) = 1;
        dataP11x4(data(:, s) == 11 & X1(:, s) == -1 & X2(:, s) == -1, s) = 1;
        dataP10x1(data(:, s) == 10 & X1(:, s) == 1 & X2(:, s) == 1, s) = 1;
        dataP10x2(data(:, s) == 10 & X1(:, s) == 1 & X2(:, s) == -1, s) = 1;
        dataP10x3(data(:, s) == 10 & X1(:, s) == -1 & X2(:, s) == 1, s) = 1;
        dataP10x4(data(:, s) == 10 & X1(:, s) == -1 & X2(:, s) == -1, s) = 1;
    end
    
    % since the covariates come from uniform discrete distribution, each
    % configuration of covariate has probability 0.25
    pX = 0.25;
    
    % Define moment (in)equalities data (n by S).
    mData_1qx1 = dataP00x1 - modelP00*pX;
    mData_1qx2 = dataP00x2 - modelP00*pX;
    mData_1qx3 = dataP00x3 - modelP00*pX;
    mData_1qx4 = dataP00x4 - modelP00*pX;
    mData_2qx1 = dataP11x1 - modelP11*pX;
    mData_2qx2 = dataP11x2 - modelP11*pX;
    mData_2qx3 = dataP11x3 - modelP11*pX;
    mData_2qx4 = dataP11x4 - modelP11*pX;
    mData_3qx1 = dataP10x1 - modelP10_lb*pX;
    mData_3qx2 = dataP10x2 - modelP10_lb*pX;
    mData_3qx3 = dataP10x3 - modelP10_lb*pX;
    mData_3qx4 = dataP10x4 - modelP10_lb*pX;
    mData_4qx1 = modelP10_ub*pX - dataP10x1;
    mData_4qx2 = modelP10_ub*pX - dataP10x2;
    mData_4qx3 = modelP10_ub*pX - dataP10x3;
    mData_4qx4 = modelP10_ub*pX - dataP10x4;
    % Stack up moments for each simulation Note: inequalities should appear first
    mData = zeros(n, 16, S);
    mbar_std = zeros(16, S);
    xi = zeros(16, S);
    epsilon = 0.000001; % introduces this parameter to avoid division by zero
    for s = 1:S
        mData(:, :, s) = [mData_3qx1(:, s), mData_3qx2(:, s), mData_3qx3(:, s), mData_3qx4(:, s),...
            mData_4qx1(:, s), mData_4qx2(:, s), mData_4qx3(:, s), mData_4qx4(:, s),...
            mData_1qx1(:, s), mData_1qx2(:, s), mData_1qx3(:, s), mData_1qx4(:, s),...
            mData_2qx1(:, s), mData_2qx2(:, s), mData_2qx3(:, s), mData_2qx4(:, s)];
        % compute studentized sample averages of mData
        mbar_std(:, s)  = mean(mData(:, :, s))./(std(mData(:, :, s)) + epsilon);
        % Additional parameter needed in DR, PR, and MR inference
        xi(:, s) = (1/kappa)*sqrt(n)*mbar_std(:, s); % Slackness measure in GMS
    end
elseif strcmp(mtype, 'agg')
    % Construct aggregate sample moments
    dataP00 = zeros(n, S); % placeholders
    dataP11 = zeros(n, S);
    dataP10 = zeros(n, S);
    for s = 1:S % loop over each simulation
        % From the data, count possible outcomes 
        % (0, 0), (1, 1) and (1, 0)
        dataP00(data(:, s) == 0, s) = 1;
        dataP11(data(:, s) == 11, s) = 1;
        dataP10(data(:, s) == 10, s) = 1;
    end
    % Define moment (in)equalities data (n by S).
    mData_1q = dataP00 - modelP00;
    mData_2q = dataP11 - modelP11;
    mData_3q = dataP10 - modelP10_lb;
    mData_4q = modelP10_ub - dataP10;
    % Stack up moments for each simulation Note: inequalities should appear first
    mData = zeros(n, 4, S); % placeholders
    mbar_std = zeros(4, S);
    xi = zeros(4, S);
    epsilon = 0.000001; % introduces this parameter to avoid division by zero
    for s = 1:S
        mData(:, :, s) = [mData_3q(:, s), mData_4q(:, s), mData_1q(:, s), mData_2q(:, s)];
        % compute studentized sample averages of mData
        mbar_std(:, s)  = mean(mData(:, :, s))./(std(mData(:, :, s)) + epsilon);
        % Additional parameter needed in DR, PR, and MR inference
        xi(:, s) = (1/kappa)*sqrt(n)*mbar_std(:, s); % Slackness measure in GMS
    end
end
end

