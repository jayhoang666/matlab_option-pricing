# matlab_option-pricing
In this respository, i will upload some Matlab code about option pricing i did, in some of which Montecarlo simulation will be used.
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
European call option under Heston model with Euler-Maruyama-Method and Milstein method

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
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Put lattice

%function price = AmputLattice(S0,K,r,T,sigma,N)
% inialization of parameters(Cox, Ross, and Rubinstein).
S0 = 10;
K = 1;
r = 0.1;
T = 4;
sigma = 1;
N = 3;
deltaT = T/N;
u = exp(sigma * sqrt(deltaT));
d = 1/u;
p = (exp(r*deltaT) - d)/(u-d);
discount = exp(-r*deltaT);
p_u = discount*p;
p_d = discount*(1-p);
% set up S values
SVals = zeros(2*N+1,1);
SVals(N+1) = S0;
for i = 1:N
SVals(N+1+i) = u*SVals(N+i);
SVals(N+1-i) = d*SVals(N+2-i);
end
% set up terminal values
PVals = zeros(2*N + 1,1);
for i =1:2:2*N+1
PVals(i) = max(K-SVals(i),0);
end
% work backwards
for tau = 1:N
for i = (tau+1):2:(2*N+1-tau)
hold = p_u*PVals(i+1) + p_d*PVals(i-1);
PVals(i) = max(hold, K-SVals(i));
end
end
price = PVals(N+1)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

European call lattice
function [price, lattice] = Ercall(S0,sigma,N,T,r,K)
deltaT = T/N;
u = exp(sigma*sqrt(deltaT));
d = 1/u;
p = ((1 + r*deltaT)-d)/(u-d);
lattice = zeros(N+1, N+1);
discount = exp(-deltaT*r);
for i = 0:N
lattice(i+1,N+1) = max(S0*(u^i)*d^(N-i)-K,0);
end

for j = N-1:-1:0
 for i = 0:j
 lattice(i+1,j+1) = discount*( p*lattice(i+2,j+2) + (1-p)*lattice(i+1,j+2));
 end
end
price= lattice(1,1)

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%Heston method for call 
function call = HestonCallQuad(kappa,theta,sigma,rho,v0,r,T,s0,K)
warning off;
call = s0*HestonP(kappa,theta,sigma,rho,v0,r,T,s0,K,1) - ...
K*exp(-r*T)*HestonP(kappa,theta,sigma,rho,v0,r,T,s0,K,2);
function ret = HestonP(kappa,theta,sigma,rho,v0,r,T,s0,K,type)
ret = 0.5 + 1/pi*quadl(@HestonPIntegrand,0,100,[],[],kappa, ...
theta,sigma,rho,v0,r,T,s0,K,type);
function ret = HestonPIntegrand(phi,kappa,theta,sigma,rho, ...
v0,r,T,s0,K,type)
ret = real(exp(-i*phi*log(K)).*Hestf(phi,kappa,theta,sigma, ...
rho,v0,r,T,s0,type)./(i*phi));
function f = Hestf(phi,kappa,theta,sigma,rho,v0,r,T,s0,type);
if type == 1
u = 0.5;
b = kappa - rho*sigma;
else
u = -0.5;
b = kappa;
end
a = kappa*theta; x = log(s0);
d = sqrt((rho*sigma*phi.*i-b).^2-sigma^2*(2*u*phi.*i-phi.^2));
g = (b-rho*sigma*phi*i + d)./(b-rho*sigma*phi*i - d);
C = r*phi.*i*T + a/sigma^2.*((b- rho*sigma*phi*i + d)*T - ...
2*log((1-g.*exp(d*T))./(1-g)));
D = (b-rho*sigma*phi*i + d)./sigma^2.*((1-exp(d*T))./ ...
(1-g.*exp(d*T)));
f = exp(C + D*v0 + i*phi*x);





