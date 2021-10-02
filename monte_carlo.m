function [] = monte_carlo(ncores,step,varargin)
% 10/20/2020 Shuowen Chen and Hiroaki Kaido
%% Inputs %%
% ncores: # of CPUs used for parallelization
% step  : 'first_stage' 'second_stage_est', 'second_stage_true', or 'power'
% Choose one of them to decide which part of the program to debug.
%
% This file conducts Monte Carlo Simulations for project "Robust Test Score
% for Incomplete Models with Nuisance Parameters"
%% Initialization
% clear
% close all
% clc

%% Setting Parameters
% beta null
beta0 = 0;
% In each sample we will draw X from this vector uniformly with replacement
covariate = [1, -1];
% true value of nuisance parameters
delta = [0.25, 0.25];

%% Random Number Generation %%
rng(123);
S = 100; % number of simulations for each sample size
parpool(ncores)


%% Standard deviations and finite-sample distribution of delta_hat

if strcmp(step,'first_stage')
    % We fix beta_alt and generate data of different sample sizes. We
    % calculate the standard deviations of beta_hat in each sample size
    % scenario (in each scenario we run MC S times).
    N = [200, 500, 1000, 2000, 5000, 10000]; % six different sample sizes
    
    % pre-allocation
    delta_BHHH = zeros(2*length(N), S);
    delta_BFGS = zeros(2*length(N), S);
    delta_LM = zeros(2*length(N), S);
    est1 = zeros(2, S);
    est2 = zeros(2, S);
    est3 = zeros(2, S);
    %parpool(8) % for parallelization (adjust this depending on # of cores)
    for j = 1:length(N) % loop over sample sizes
        n = N(j);
        % Call data_generation function to generate data
        [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, ...
            delta, n, S, 'Complete');
        % loop over simulations
        parfor i = 1:S
            % Conduct Restricted MLE
            est1(:,i) = rmle([0,0], n, 'BHHH', data_com(:,i), x1_com(:,i), x2_com(:,i));
            est2(:,i) = rmle([0,0], n, 'BFGS', data_com(:,i), x1_com(:,i), x2_com(:,i));
            est3(:,i) = rmle([0,0], n, 'LM-BFGS', data_com(:,i), x1_com(:,i), x2_com(:,i));
        end
        % Store estimates
        delta_BHHH((2*j-1):(2*j), :) = est1;
        delta_BFGS((2*j-1):(2*j), :) = est2;
        delta_LM((2*j-1):(2*j), :) = est3;
    end
    
    %% Checking the first stage results (BHHH, BFGS, L-BFGS-B)
    % count number of NAs across sample sizes
    missing_BHHH = sum(isnan(delta_BHHH),2);
    missing_BFGS = sum(isnan(delta_BFGS),2);
    missing_LM = sum(isnan(delta_LM),2);
    % Checking how close the distribution looks like normal
    normplot(delta_BHHH(2,:)) % small sample is problematic, need to revise
    normplot(delta_BHHH(11,:))
    normplot(delta_BHHH(12,:)) % the tail is not satisfying
    J = length(N);
    
    for k = 1:3
        figure(k)
        for j=1:J
            subplot(J,2,2*j-1)
            switch k
                case 1
                    histogram(delta_BHHH(2*(j-1)+1,:),floor(2*S^(1/3)));
                case 2
                    histogram(delta_BFGS(2*(j-1)+1,:),floor(2*S^(1/3)));
                case 3
                    histogram(delta_LM(2*(j-1)+1,:),floor(2*S^(1/3)));
            end
            title(['n=' num2str(N(j))])
            subplot(J,2,2*j)
            switch k
                case 1
                    histogram(delta_BHHH(2*j,:),floor(2*S^(1/3)));
                case 2
                    histogram(delta_BFGS(2*j,:),floor(2*S^(1/3)));
                case 3
                    histogram(delta_LM(2*j,:),floor(2*S^(1/3)));
            end
            title(['n=' num2str(N(j))])
        end
        switch k
            case 1
                sgtitle('Distribution of RMLE (BHHH)')
                saveas(gcf,'../Figures/histogram_BHHH.pdf')
            case 2
                sgtitle('Distribution of RMLE (BFGS)')
                saveas(gcf,'../Figures/histogram_BFGS.pdf')
            case 3
                sgtitle('Distribution of RMLE (L-BFGS-B)')
                saveas(gcf,'../Figures/histogram_LBFGSB.pdf')
        end
    end
    
    %% Size properties
elseif strcmp(step,'second_stage_est')
    % generate data under the null
    level = 0.05; % nominal size
    lambda = 0.01; % regularization parameter
    % critical value using the limiting distribution
    %     [cv,simstat] = crt(level);
    N = [2500, 5000, 7500];
    %     N=4000;
    J = length(N);
    X = [1,1; 1,-1; -1,1; -1,-1]; % combination of covariates
    % pre-allocation to store g_delta statstics and sup test statistic
    gn = zeros(2*length(N), S);
    test = zeros(length(N), S);
    outcome_com = zeros(16,S,length(N));
    delta_hat_com1 = zeros(length(N),S);
    delta_hat_com2 = zeros(length(N),S);
    cv = zeros(length(N),S);
    for j = 1:length(N) % loop over sample size
        n = N(j);
        % Call data_generation function to generate data
        [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, ...
            delta, n, S, 'Complete');
        outcome_com(:,:,j) = counting(data_com, x1_com, x2_com);
        % loop over number of simulations
        gn_temp = zeros(2,S);
        parfor i = 1:S
            % Conduct Restricted MLE
            delta_hat_com = rmle([0,0], n, 'BFGS', data_com(:,i), x1_com(:,i), x2_com(:,i));
            delta_hat_com1(j,i) = delta_hat_com(1);
            delta_hat_com2(j,i) = delta_hat_com(2);
            % generate var-cov matrices and gn statistics (generated under the null beta=0)
            %[test(j,i),gn_temp(:,i), varI] = stat(delta_hat_com, X, outcome_com(:,i,j), n,lambda);
            [test(j,i),gn_temp(:,i), varI] = stat_ab(delta_hat_com, X, outcome_com(:,i,j), n);
            cv(j,i) = crt_ab(level,varI);
        end
        gn((2*j-1):(2*j),:) = gn_temp;
        filename = ['../Results/Matfiles/test_crt_n' num2str(n) '_S' num2str(S) '.mat'];
        save(filename)
    end
    % count missing gn values, can see coincides with missings from delta_hat
    gn1 = gn(1:2:2*J,:);
    gn2 = gn(2:2:2*J,:);
    missing2 = sum(isnan(gn1),2);
    
    
    % plot gn statistics
    figure(1)
    for j = 1:J
        subplot(J,2,2*j-1)
        histogram(gn1(j,:),floor(2*S^(1/3)));
        title(['n=' num2str(N(j))])
        subplot(J,2,2*j)
        histogram(gn2(j,:),floor(2*S^(1/3)));
        title(['n=' num2str(N(j))])
    end
    %sgtitle('Distribution of gn under H0')
    %saveas(gcf,'../Figures/dist_gn_est.pdf')
    
    figure(2)
    for j = 1:J
        subplot(J,2,2*j-1)
        normplot(gn1(j,:));
        title(['n=' num2str(N(j))])
        subplot(J,2,2*j)
        normplot(gn2(j,:));
        title(['n=' num2str(N(j))])
    end
    %sgtitle('Probability plot of gn under H0')
    %saveas(gcf,'../Figures/pplot_gn_est.pdf')
    
    figure(3)
    for j=1:J
        subplot(J,1,j)
        histogram(test(j,:))
        title(['n=' num2str(N(j))])
    end
    sgtitle('Distribution of the test statistic under H0')
    saveas(gcf,'../Figures/dist_stat_est.pdf')
    
    % figure(4)
    % for j=1:J
    %     subplot(J,1,j)
    %     qqplot(sim_stat,test(j,:));
    %     title(['n=' num2str(N(j))])
    % end
    %sgtitle('Probability plot of the test statistic under H0')
    %saveas(gcf,'../Figures/pplot_stat_est.pdf')
    
    % compute the real size
    [real, l] = actualsize(test, N, cv);
    % From l2, we see that without regularization there are many missing value
    % for smalle samples
    
    %% The previous simulation uses delta_hat, now use true delta value
elseif strcmp(step,'second_stage_true')
    % generate data under the null
    level = 0.05; % nominal size
    lambda = 0.1; % regularization parameter
    % critical value using the limiting distribution
    [cv,simstat] = crt(level);
    N = 20000;
    J = length(N);
    X = [1,1; 1,-1; -1,1; -1,-1]; % combination of covariates
    % pre-allocation to store g_delta statstics and sup test statistic
    gn = zeros(2*length(N), S);
    test = zeros(length(N), S);
    outcome_com = zeros(16,S,length(N));
    for j = 1:length(N) % loop over sample size
        n = N(j);
        % Call data_generation function to generate data
        [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, ...
            delta, n, S, 'Complete');
        outcome_com(:,:,j) = counting(data_com, x1_com, x2_com);
        % loop over number of simulations
        gn_temp = zeros(2,S);
        parfor i = 1:S
            % generate var-cov matrices and gn statistics (generated under the null beta=0)
            [test(j,i),gn_temp(:,i)] = stat(delta', X, outcome_com(:,i,j), n,lambda);
        end
        gn((2*j-1):(2*j),:) = gn_temp;
    end
    % count missing gn values, can see coincides with missings from delta_hat
    gn1 = gn(1:2:2*J,:);
    gn2 = gn(2:2:2*J,:);
    missing2 = sum(isnan(gn1),2);
    % plot gn statistics
    figure(1)
    for j = 1:J
        subplot(J,2,2*j-1)
        histogram(gn1(j,:),floor(2*S^(1/3)));
        title(['n=' num2str(N(j))])
        subplot(J,2,2*j)
        histogram(gn2(j,:),floor(2*S^(1/3)));
        title(['n=' num2str(N(j))])
    end
    %sgtitle('Distribution of gn under H0')
    %saveas(gcf,'../Figures/dist_gn_true.pdf')
    
    figure(2)
    for j = 1:J
        subplot(J,2,2*j-1)
        normplot(gn1(j,:));
        title(['n=' num2str(N(j))])
        subplot(J,2,2*j)
        normplot(gn2(j,:));
        title(['n=' num2str(N(j))])
    end
    %sgtitle('Probability plot of gn under H0')
    %saveas(gcf,'../Figures/pplot_gn_true.pdf')
    
    figure(3)
    for j=1:J
        subplot(J,1,j)
        histogram(test(j,:))
        title(['n=' num2str(N(j))])
    end
    %sgtitle('Distribution of the test statistic under H0')
    %saveas(gcf,'../Figures/dist_stat_true.pdf')
    
    figure(4)
    for j=1:J
        subplot(J,1,j)
        qqplot(simstat,test(j,:));
        title(['n=' num2str(N(j))])
    end
    %sgtitle('Probability plot of the test statistic under H0')
    %saveas(gcf,'../Figures/pplot_stat_true.pdf')
    
    % compute the real size
    [real, l] = actualsize(test, N, cv);
    
    
    %% Power analysis
elseif strcmp(step,'power')
    %%
    level = 0.05; % nominal size
    n = 7500;
    X = [1,1; 1,-1; -1,1; -1,-1]; % combination of covariates
    DGP = 'LFP';
    
    % pre-define alternatives
    h_alt  = -(eps:0.5:15)';
    h_alt2 = -(eps:0.5:15)';
    K = length(h_alt);
    % pre-allocation to store sup test statistic for each selection mechanism
    test = zeros(K, S);
    cv   = zeros(K,S);
    for k=1:K
        % Define true beta that generates the data (Change between beta_alt
        % and beta_alt2 in the data_generation function)
        beta_alt = h_alt(k,:)/sqrt(n);
        beta_alt2 = h_alt2(k,:)/sqrt(n);
        % Call data_generation function to generate data
        [data, x1, x2] = data_generation([beta_alt,beta_alt2], covariate, ...
            delta, n, S, DGP);
        outcome = counting(data, x1, x2);
        % loop over number of simulations
        parfor i = 1:S
            % Conduct Restricted MLE
            delta_hat = rmle([0,0], n, 'BFGS', data(:,i), x1(:,i), x2(:,i));
            % generate var-cov matrices and gn statistics (generated under the null beta=0)
            [test(k,i), ~,varI_iid] = stat_ab(delta_hat, X, outcome(:,i), n);
            cv(k,i) = crt_ab(level,varI_iid);
        end
    end
    filename = ['../Results/Matfiles/test_power_DGP' DGP '_n' num2str(n) '_S' num2str(S) '.mat'];
    save(filename)
    
    % Post Model Selection Inference
elseif strcmp(step,'ci')
    k = varargin{1};
    addpath(genpath('./KMS'))
    % This part presumes you ran the 'power' step already.
    level = 0.05; % nominal size
    n = 2500;
    X = [1,1; 1,-1; -1,1; -1,-1]; % combination of covariates
    DGP = 'IID';
    comp = 1; % component of interest
    
    % set shrinkage factor
    kappa = 1/sqrt(log(n));
    
    % pre-define alternatives
    betaval = [-0.2,-0.1,-0.05,-0.025,-eps];
    % pre-allocation to store sup test statistic for each selection mechanism
    test = zeros(S,1);
    cv   = zeros(S,1);
    citype = zeros(S,1); % =1 (KMS-ci), =0 (Wald-ci)
    ciu  = zeros(S,1);
    cil  = zeros(S,1);
    output = cell(S,1);
    % True beta
    beta_alt = betaval(k);
    beta_alt2 = betaval(k);
    % Call data_generation function to generate data
    [data, x1, x2] = data_generation([beta_alt,beta_alt2], covariate, ...
        delta, n, S, DGP);
    outcome = counting(data, x1, x2);
    % loop over number of simulations
    for i = 1:S
        % Conduct Restricted MLE
        delta_hat = rmle([0,0], n, 'BFGS', data(:,i), x1(:,i), x2(:,i));
        % generate var-cov matrices and gn statistics (generated under the null beta=0)
        [test(i), ~,varI_iid] = stat_ab(delta_hat, X, outcome(:,i), n);
        cv(i) = crt_ab(level,varI_iid);
        
        if test(i) > kappa*cv(i)
            citype(i) = 1;
            [cil(i),ciu(i),output{i}] = get_KMS(data(:,i),x1(:,i),x2(:,i),level,comp);
        else
            citype(i) = 0;
            [cil(i),ciu(i)] = get_wald(delta_hat,X,outcome(:,i),n,level,comp);
        end
    end
    filename = ['../Results/Matfiles/ci_DGP' DGP '_k' num2str(k) '_n' num2str(n) '_S' num2str(S) '.mat'];
    save(filename)
    
end
%date = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
%filename = strcat('../Results/Results_',step,date,'.mat');
%save(filename)
end