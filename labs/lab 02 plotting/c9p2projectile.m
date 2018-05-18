%c9p2projectile.m
%This script produces a 3D graph of a projectile fired into the air given
%the initial boundary conditions of position, velocity, and acceleration.
%East is in the +x direction, and North is in the +y direction. All units
%are SI units.

g=[0,0,9.81];
r0=[0,0,0];

v0=250;%input('Enter initial speed       v0 (m/s): ');
alt=pi/180*65;%input('Enter the initial altitude angle relative to horizon (deg): ');
az=pi/180*0;%input('Enter the initial azimuth angle relative to North (deg): ');
interval=10;%input('Enter the interval of motion (s): ');
increment=.1;%input('Enter the time increment     (s): ');
windV=[-65,0,0];%input('Enter the wind velocity [vx,vy,vz]: ');

v0=[v0*sin(az),v0*cos(az),v0*cos(alt)]+windV;

t=0:increment:interval;
for I=1:length(t),
    R(1,I)=r0+v0*t(1,I)+.5*g*t(1,I)^2;
end
disp(R);
disp(t);
plot(t,R);
