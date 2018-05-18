function []=de2nd(dt,tmin,tmax)

%tmin=0;
%tmax=5;
%dt=0.01;
t=tmin;
yi=1;
zi=2;
zpi=0;

ts=tmin;
ys=yi;

while t <= tmax
    yi=y(dt,yi,zi);
    zi=z(dt,zi,zpi);
    zpi=zp(t,yi,zi);
    
    t=t+dt/2;
    ts=[ts t];
    ys=[ys yi];
end


plot(ts,ys);
title('Modified Euler 2nd Order');
xlabel('t values');
ylabel('y values');
end

function v=y(h,yi,ypi)
v=yi+h/2*ypi;
end

function v=z(h,zi,zpi)
v=zi+h/2*zpi;
end

function v=zp(t,y,z)
v=-2*z-y+exp(-t);
end