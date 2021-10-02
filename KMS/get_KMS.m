function [cil, ciu, KMS_output] = get_KMS(y,xx1,xx2,level,comp)
%% Key Parameters:
method      = 'KMS';    % Method - either AS or KMS
DGP         = 11;        % DGP9 is for the empirical example
alpha       = level;     % Significance level
component   = comp;        % Component of theta to build confidence interval around
phi = NaN;
% phi   = @(x)(min(x,0));

%% Load data
KMSoptions_app  = KMSoptions();
suppX = [1 -1 1 -1 ; 1 -1 1 1 ; 1 1 1 -1 ; 1 1 1 1];

[n,~] = size(y);
y1 = 1*(y>=10);
y2 = 1*((y==1)+(y==11));
W = [y1, y2, ones(n,1), xx1, ones(n,1), xx2];

dim_suppX = size(suppX,1);
psuppX = ones(dim_suppX,1)./dim_suppX;
% psuppX = zeros(dim_suppX,1);
% for i=1:n
%     for j = 1:dim_suppX
%         if W(i,3:6)==suppX(j,:)
%             psuppX(j) = psuppX(j)+1;
%         end
%     end
% end
% psuppX = psuppX./n;

KMSoptions_app.suppX = suppX;
KMSoptions_app.psuppX = psuppX;

% Set other KMSoptions_app
KMSoptions_app.DGP = DGP;
KMSoptions_app.n = n;
KMSoptions_app.component = component;
KMSoptions_app.BCS_EAM = 0;
KMSoptions_app.boundary = 1;
%KMSoptions_app.boundary = 0;
seed = KMSoptions_app.seed;
B    = KMSoptions_app.B;
stream = RandStream('mlfg6331_64','Seed',seed);
RandStream.setGlobalStream(stream)

%% Parameters
type = 'two-sided';         % Two-sided or one sided test?  Set to 'one-sided-UB' or 'one-sided-LB' or 'two-sided'
kappa =NaN;                 % Default kappa function

dim_p=4;
LB_theta = [-2;-2;-4;-4];  % Lower bound on parameter space
UB_theta = [2;2;0;0];
theta_0 = 0.5*UB_theta + 0.5*LB_theta;
p = zeros(dim_p,1);
p(component) = 1;
KMSoptions_app.S =  0;   % Rho Polytope Constraints
A_theta = [];
b_theta = [];
CVXGEN_name = strcat('csolve_DGP',num2str(DGP));

%% Run KMS
[KMS_confidence_interval,KMS_output] = KMS_0_Main(W,theta_0,...
            p,[],LB_theta,UB_theta,A_theta,b_theta,alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions_app);
cil = KMS_confidence_interval(1);
ciu = KMS_confidence_interval(2);
end
