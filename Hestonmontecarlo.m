% European call option under Heston model with Euler-Maruyama-Method and Milstein method

clear all, randn('state',100)
M = 50000;                        % number of paths
N = 100;                        % number of time steps
T = 1; % maturity
h = T/N;   %delta t
S0 = 100; sigma20 = 0.0625; K=100; % initial values
kappa = 2; theta = 0.4; nu = 0.2; rho = -0.7;
r=0.02;

% two dimensional Brownian motion
dW1 = randn(M,N+1)*sqrt(h);
dW2 = rho*dW1 + sqrt(1-rho^2)*randn(M,N+1)*sqrt(h);

% Initialisation of stock price and volatility for Euler Maruyama method
S = S0*ones(M,N+1);
sigma2 = sigma20*ones(M,N+1);

% Solution of the SDE-System with Euler-Maruyama-Method
for i = 1:N
    sigma2(:,i+1) = sigma2(:,i) + kappa*(theta-sigma2(:,i))*h ...
        + nu*sqrt(sigma2(:,i)).*dW2(:,i);
    S(:,i+1) = S(:,i).*(1 + r*h + sqrt(sigma2(:,i)).*dW1(:,i));
end

payoff = max(0,S(:,end)-K);
V = exp(-r*T)*mean(payoff)
% Initialisation of stock price and volatility for Milstein method
S1 = S0*ones(M,N+1);
sigma3 = sigma20*ones(M,N+1);
% Solution of the Milstein method

for i = 1:N
    sigma3(:,i+1) = sigma3(:,i) + kappa*(theta-sigma3(:,i))*h ...
        + nu*sqrt(sigma3(:,i)).*dW2(:,i)+1/4*nu^2.*(dW2(:,i).^2-h);
    S1(:,i+1) = S1(:,i).*(1 + r*h + sqrt(sigma3(:,i)).*dW1(:,i))+1/2*sigma3(:,i).*S1(:,i).*(dW1(:,i).^2-h);
end

payoff1 = max(0,S1(:,end)-K);
V1 = exp(-r*T)*mean(payoff1)