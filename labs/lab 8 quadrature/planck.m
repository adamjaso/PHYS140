%planck.m
%This script determines the value visible fraction of useful energy from
%the planck blackbody equation
function []=planck()
%define constants
k=1.38e-23;
T=2000;%Kelvin
h=6.625e-34;
c=2.9979e8;
xToL=h*c/(k*T);

%integration constant
A=(2*pi*k^4*T^4)/(h^3*c^2);

%table values
xtab=[-0.99313 -0.96397 -0.91223 -0.83912 -0.74633 -0.63605 -0.51087 -0.37371 -0.22779 -0.07653 .07653 .22779 .37371 .51087 .63605 .74633 .83912 .91223 .96397 .99313];
wtab=[0.01761 0.04060 0.06267 0.08328 0.10193 0.11819 0.13169 0.14210 0.14917 0.15275 .15275 .14917 .14210 .13169 .11819 .10193 .08328 .06267 .04060 .01761];

Etop=Enum(A,xToL,xtab,wtab);
Ebot=Etot(A,xtab,wtab);

fprintf('Numerator   = %f J\n',Etop);
fprintf('-----------------------------------\n');
fprintf('Denominator = %f J\n',Ebot);
fprintf('The Ratio = %f\n',Etop/Ebot);
end

%calculates the total energy
function eval=Etot(A,xtab,wtab)
Ef=Efirst(xtab,wtab);
Es=Esecond(xtab,wtab);
eval=-A*(Ef+Es);
fprintf('Etot: first  integral = %f\n',Ef);
fprintf('Etot: second integral = %f\n',Es);
end
%calculates the finite part of the total energy integral
function I=Efirst(xtab,wtab)

a=0;
b=3;
s1=(b-a)/2;
s0=(b+a)/2;
Wi=(b-a)/2*wtab;
zi=s1*xtab+s0;

I=0;
for i=1:length(zi)
    I=I+Wi(i)*f1(zi(i));
end

end
%calculates the inifinite part of the total energy integral
function I=Esecond(xtab,wtab)

a=1/3;
b=0;
s1=(b-a)/2;
s0=(b+a)/2;
Wi=(b-a)/2*wtab;
zi=s1*xtab+s0;

I=0;
for i=1:length(zi)
    I=I+Wi(i)*f2(zi(i));
end
I=-I;
end

%integrand function
function eval=f1(x)
eval=x^3/(exp(-x)-1);
end

%Esecond integrand
function eval=f2(x)
eval=x^-5*(exp(1/x)-1)^-1;
end

%calculates the numerator energy
function I=Enum(A,xToL,xtab,wtab)
a=xToL*1/400e-9;
b=xToL*1/700e-9;
s1=(b-a)/2;
s0=(b+a)/2;
zi=s1*xtab+s0;
Wi=(b-a)/2*wtab;

I=0;
for i=1:length(zi)
    I=I+Wi(i)*f1(zi(i));
end

I=A*I;
end

