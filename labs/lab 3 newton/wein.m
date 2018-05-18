%wein.m
%This program determines the wavelength that maximizes Planck's black body
%equation

h=6.6262e-34;
c=3e8;
k=1.3806e-23;
x=10;
dx = newton(x);
err=0.0000000001;

while abs(dx) > err
    dx = newton(x);
    x = x + dx;
end

constant=(h*c)/(x*k);
fprintf('x = %2.10f\n',constant);

