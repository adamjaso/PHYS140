%trap52.m
%This script determines the mass of the sun from a model using the trapezoidal rule

%1.98892 × 10^30

Rsun = 6.955e8;
rs = 0:.1:1;
rhos = [160000 90000 40000 13000 4000 1000 400 80 20 2 .0003];
M = 4 * pi * Rsun^2 * sum((rs.^2 * .1) .* rhos);

r2xrhos = rs.^2 .* rhos;

I=0;
for i=2:length(rs)
    I = I + .5 * ((rs(i)-rs(i-1))*(r2xrhos(i)+r2xrhos(i-1)));
end

fprintf('M = %e\n',4*pi*Rsun^3*I);
