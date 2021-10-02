function [] = Simulation_BCS(ncores,part)
% 5/10/2021 Shuowen Chen and Hiroaki Kaido
% Monte Carlo for comparison with Bugni, Canay and Shi (2017)
%% 0. Parameters
% seed for replication
rng(123);
parpool(ncores)
% In each sample we will draw X from this vector uniformly with replacement
covariate = [1, -1];
% true value of nuisance parameters
delta = [2, 1.5];
% number of simulations
S = 1999;
% N = [200, 500, 1000, 2000, 5000, 10000]; % sample sizes
n = 7500;
% # of bootstrap draws
B = 500;
% shocks \xi in eq (2.8) pp. 5
W2_AA = norminv(rand(B,n));
% constant in the paper;
C_list = [1];
% GMS Thresholding parameter, preferred choice Eq. 4.4 in Andrews and Soares (2010).
kappa_list = C_list*sqrt(log(n));
% %% 1. DGP under the null
% [data_com, x1_com, x2_com] = data_generation([0, 0], covariate, delta, n, S, 'Complete');
%% 2. Implement Test MR in BCS (algorithm 2.1)
options = optimset('Display','off','Algorithm','interior-point'); % from BCS
% significance level
DGP = 'LFP';
alpha = 0.1;
% Placeholders
minQn_DR = NaN(size(kappa_list,1), B);
minQn_PR = NaN(size(kappa_list,1), B);
minQn_MR = NaN(size(kappa_list,1), B);
cn_MR  = NaN(1,size(kappa_list,1));

h_alt  = -(eps:0.5:15)';
h_alt2 = -(eps:0.5:15)';
K = length(h_alt);
Tn_MRsim = NaN(K,S);
cn_MRsim = NaN(K,S);
pbegin = (part-1)*4+1;
pend = part*4;
if part == 8
    pend = 31;
end
% Implementation starts here
for k = pbegin:pend
    beta_alt = h_alt(k,:)/sqrt(n);
    beta_alt2 = h_alt2(k,:)/sqrt(n);
    % Call data_generation function to generate data
    [data_com, x1_com, x2_com] = data_generation([beta_alt,beta_alt2], covariate, ...
        delta, n, S, DGP);
    for s = 1:S % loop over simulations
        disp(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the test statstic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Pick a large value as initial function value for the minimization
        min_value = 10^10; % from BCS
        % Multiple starting points using Latin Hypercube
        lb = delta - 1;
        ub = delta + 1;
        M = 49; % number of starting points using lhsscale
        start = [delta; lhsscale(M, 2, ub, lb)]; % 1st starting point is oracle
        % placeholder for minimization results
        min_outcomes = zeros(M+1, 3);
        
        % loop over starting points for minimization
        for initial = 1:(M + 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (Line 20 in BCS Algorithm 2.1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Instead of set minimizer we have point estimator of nuisance par
            try
                [delta_est,Qn_aux,bandera] =  fminunc(@(x) Qn_function_X(x, [0, 0], data_com(:,s), x1_com(:,s), x2_com(:,s), 'agg'), start(initial, :), options);
                % check whether minimization is successful and reduced value
                if Qn_aux  < min_value && bandera >= 1
                    Qn_minimizer = delta_est;
                    min_value = Qn_aux;
                end
            catch
                min_value = NaN;
                bandera = 0;
            end
            
            % if minimization is successful, collect minimizer and its value;
            if bandera >= 1
                min_outcomes(initial,:) = [min_value, Qn_minimizer];
            end
        end
        % return the test statstic
        minQn = min_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Test MR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kappa_type =1;
        parfor b = 1:B % loop over bootstrap draws
            % (1) DR test (equivalent to plugging in nuisance estimator and null value)
            % compute simulated DR criterion function
            minQn_DR(kappa_type, b) = Qn_MR_function_X(delta_est, [0,0], data_com(:,s), x1_com(:,s), x2_com(:,s), kappa_list(kappa_type), W2_AA(b,:), 1, 'agg');
            
            % (2) PR test
            [delta_aux,min_value_PR] =  fmincon(@(x) Qn_MR_function_X(x, [0,0], data_com(:,s), x1_com(:,s), x2_com(:,s), kappa_list(kappa_type), W2_AA(b,:), 2, 'agg'),...
                delta,[],[],[],[],...
                [1,1],[2.5,5],[],options);
            % compute simulated PR criterion function
            minQn_PR(kappa_type, b) = min_value_PR ;
            
            % (3) MR test
            minQn_MR(kappa_type, b) = min(minQn_DR(kappa_type, b), minQn_PR(kappa_type, b));
        end
        % Compute critical values
        cn_MR(kappa_type) = quantile(minQn_MR(kappa_type,:), 1-alpha);
        cn_MRsim(k,s) = cn_MR;
        Tn_MRsim(k,s) = minQn;
        %     accept_testMR = repmat(minQn,1,size(kappa_list,1)) <= cn_MR;
    end
end
filename = ['../Results/Matfiles/BCS_power_DGP' DGP '_n' num2str(n) '_S' num2str(S) 'part' num2str(part) '.mat'];
save(filename)

end
%%
% % placeholder
% delta_hat = NaN(2, S);
% for s = 1:S
%     delta_hat(:, s) = fminunc(@(x) Qn_function_X(x, [0, 0], data_com(:,s), x1_com(:,s), x2_com(:,s), 'agg'), delta);
% end
%
% %%
% X1 = x1_com(:,1);
% X2 = x2_com(:,1);
% beta = [0, 0];
% delta0 = [2, 1.5];
% % create a grid centered around delta0
% [A, B] = meshgrid((delta0(1)-2):0.05:(delta0(1)+2), (delta0(2)-2):0.05:(delta0(2)+2));
% val = NaN(size(A));
%
% % loop over grid points to plot the objective function
% for i = 1:size(A, 1)
%     for j = 1:size(A, 2)
%         val(i,j) = Qn_function_X([A(1, i), B(j, 1)], [0, 0], data_com(:, 1), X1, X2, 'agg');
%     end
% end
% figure;
% surf(A, B, val)
% colorbar
