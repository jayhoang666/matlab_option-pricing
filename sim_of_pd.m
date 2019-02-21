
   strikes = linspace(.8,1.2,11);
 mats = linspace(.3,3,11);      %maturities
 for i = 1:11
    for j = 1:11
        price = HestonCallQuad(2,.04,.1,0.5,.04,.01,mats(i),1,strikes(j));
        prices(i,j) = price;
        Volatility(i,j) = blsimpv(1, strikes(j), 0.01 , ...
                                                mats(i), price);
    end
 end

  [strike mat] = meshgrid(strikes,mats);
 surf(mat,strike,Volatility),xlabel('Maturity(years)'), ...
 ylabel('Strike'),Title('\rho = 0.5'),zlabel('Implied Volatility');
 figure;
for i = 1:11
   for j = 1:11
       price = HestonCallQuad(2,.04,.1,0,.04,.01,mats(i),1,strikes(j));
       prices(i,j) = price;
       Volatility(i,j) = blsimpv(1, strikes(j), 0.01 , ...
                                            mats(i), price);
   end
end
surf(mat,strike,Volatility),xlabel('Maturity(years)'),
ylabel('Strike'),Title('\rho = 0'),zlabel('Implied Volatility');
figure;
for i = 1:11
   for j = 1:11
price = HestonCallQuad(2,.04,.1,-0.5,.04,.01,mats(i),...
                                          1,strikes(j));
   end
end
