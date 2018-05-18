function []=runge4(dt)
%runge4.m
%This function finds the numerical solution to the ODE: y' = 3*e^(-4*t)-2*y

tmin=0;
tmax=5;
%dt=0.01;
t=0;
ti=0;
yi=1;
k1=dt*diffeq(ti,yi);
ts=ti;
ys=yi;

while t < tmax
    t=t+dt;
    ti=t;
    k1=dt*diffeq(ti,yi);
    k2=dt*diffeq(ti+dt/2,yi+k1/2);
    k3=dt*diffeq(ti+dt/2,yi+k2/2);
    k4=dt*diffeq(ti+dt,yi+k3);
    yi=yi+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
    
    ts=[ts ti];
    ys=[ys yi];
end

plot(ts,ys);
title('Runge / Cutta: 4 - point');
xlabel('t values');
ylabel('y values');
end

function diffeq=diffeq(t,y)
diffeq=3*exp(-4*t)-2*y;
end