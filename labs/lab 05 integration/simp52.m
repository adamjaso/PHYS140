%simp52.m
%This script determines the mass of the sun from a model using Simpson's rule

Rsun = 6.955e8;
rs = 0:.1:1;
rhos = [160000 90000 40000 13000 4000 1000 400 80 20 2 .0003];

I = 0;

for i=1:(length(rs)-2)/2
    I = I + 1/3 * (rs(2*i+1) - rs(2*i)) * (r2xrho(rs(2*i),rhos(2*i)) + 4 * r2xrho(rs(2*i+1),rhos(2*i+1)) + r2xrho(rs(2*i+2),rhos(2*i+2)));
end

fprintf('M = %e\n',4*pi*Rsun^3*I);