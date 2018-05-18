%act3prob3.m
%This is a solution to Activity #3 Problem #3

xAxisLabel='Time (s)';
t=0:0.1:8;
x=-.1*t.^4+.8*t.^3+10*t-70;
v=-.4*t.^3+2.4*t.^2+10;
a=-1.2*t.^2+4.8*t;

subplot(3,1,1),plot(t,x);
xlabel(xAxisLabel);
ylabel('Distance (m)');
subplot(3,1,2),plot(t,v);
xlabel(xAxisLabel);
ylabel('Speed (m/s)');
subplot(3,1,3),plot(t,a);
xlabel(xAxisLabel);
ylabel('Acceleration (m/s^2)');