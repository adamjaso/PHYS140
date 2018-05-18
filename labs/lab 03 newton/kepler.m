%kepler.m
%This program finds the eccentric anomaly, E, of an elliptical orbit.

P=75.99;
a=17.94;
angSp=2*pi/P;
e=0.967;
t=15;
err=1e-7;

E=angSp*t+e*sin(angSp*t)+0.5*e^2*sin(2*angSp*t);
dE=knewton(E,e,angSp,t);

while abs(dE) > err
    dE=knewton(E,e,angSp,t);
    E=E+dE;
end

r=a*(1-e*cos(E));
th=180/pi*acos((cos(E)-e)/(1-e*cos(E)));
fprintf('E(t = %2.1f s) = %5.8f rad\nr(theta = %2.4f deg) = %5.8f AU\n',t,E,th,r);
