function []=eulermod(dt)
%euler.m
%This function finds the numerical solution to the ODE: y' = 3*e^(-4*t)-2*y

tmin=0;
tmax=5;
%dt=0.0005;
t=tmin;
yi=1;
ypi=diffeq(tmin,yi);

ts=tmin;
ys=yi;
while t <= tmax
    yi=apx(yi,ypi,dt);
    ypi=diffeq(t,yi);
    t=t+dt/2;
    ts=[ts t];
    ys=[ys yi];
end

plot(ts,ys);
title('Modified Euler');
xlabel('t values');
ylabel('y values');
end

function diffeq=diffeq(t,y)
diffeq=3*exp(-4*t)-2*y;
end

function apx=apx(y,yp,d)
apx=y+d/2*yp;
end