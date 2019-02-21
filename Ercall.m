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
lattice
