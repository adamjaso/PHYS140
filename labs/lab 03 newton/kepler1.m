%kepler.m
%This program finds the eccentric anomaly, E, of an elliptical orbit.

yrToSec=365*24*3600;

P=75.99;
a=17.94;
angSp=2*pi/P;
e=0.967;
t=0:1:P;
err=1e-7;
E(length(t))=0;


for n=1:length(t)
    E(n)=angSp*t(n)+e*sin(angSp*t(n))+0.5*e^2*sin(2*angSp*t(n));
    dE=knewton(E(n),e,angSp,t(n));
    c=0;
    while abs(dE) > err
        before = E(n);
        dE=knewton(E(n),e,angSp,t(n));
        E(n)=E(n)+dE;
        after = E(n);
        if before < after
            c = c + 1;
        end
        if c > 10
            break;
        end
    end
    %fprintf('E(t = %2.3f) = %5.8f\n',t(n),E(n));
end
plot(t,E);
