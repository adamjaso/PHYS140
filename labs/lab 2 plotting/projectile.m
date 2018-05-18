%projectile.m
%This script generates the plot of the path of a fired projectile moving in
%a frictionless environment above the surface of earth

alt=pi/180*65;
az=pi/180*90;
windv=[-30,0,0];
v0mag=250;
r0=[3000,0,0];
v0=[v0mag*cos(az),v0mag*sin(az),v0mag*sin(alt)]+windv;
g=[0,0,9.81];

t=0:0.1:45;
x=r0(1)+v0(1)*t-.5*g(1)*t.^2;
y=r0(2)+v0(2)*t-.5*g(2)*t.^2;
z=r0(3)+v0(3)*t-.5*g(3)*t.^2;
plot3(x,y,z);
title('The motion of a projectile in 45 s');
xlabel('+x East (m)');
ylabel('+y North (m)');
zlabel('+z Up (m)');
grid on;


