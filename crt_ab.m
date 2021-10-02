% 3/10/2020 Hiroaki Kaido and Shuowen Chen
% The function calculates the critical value of the tn statistic
% Recall that the limiting distribution of tn is the sup of bivariate
% standard normal distribution
function crit = crt_ab(nominal,varI)
rng(pi); % set seed for random number generations
J = 2;    % # of components in the score
R = 5000; % # of draws
z = randn(J, R); % J by R vector of draws from standard normal
[LI,DI] = ldl(varI);
ztilde=LI*DI^(0.5)*z;
%     DI(DI<tol) = tol; % truncate negative values

first = sum(z.^2); % limiting distribution
second = zeros(1,R);
for r=1:R
    Vmin = @(x) quadform(varI,ztilde(:,r)-x);
    options = optimoptions('fmincon','Display','off');
    [~,Vstar] = fmincon(Vmin,[0;0],eye(2),zeros(2,1),[],[],[],[],[],options);
    second(r) = Vstar;
end
sim_stat = first-second;
crit = quantile(sim_stat, (1-nominal));
end

function t = quadform(A,g)
gtilde = A^(-0.5)*g;
t = gtilde'*gtilde;
end