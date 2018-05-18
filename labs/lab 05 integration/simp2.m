%simp.m
%This script determines the interval of sin(x) using Simpson's rule

a = 0; %input('Enter the interval start: ');
b = pi/2; %input('Enter the interval end:   ');
n = 1200; %input('Enter the intervals:      ');
err = 1e-6;

xs = linspace(a,b,n);
I = 0;

while abs(I - (-cos(b)+cos(a))) > err
    I=0;
    for i=1:(n-2)/2
        I = I + 1/3 * (xs(2*i+1) - xs(2*i)) * (sin(xs(2*i)) + 4 * sin(xs(2*i+1)) + sin(xs(2*i+2)));
    end
    n = n + 100;
    xs = linspace(a,b,n);
end

if abs(I-1) < err
    fprintf('INT(sin(x),x,%f,%f) = %1.10f\nerr < %1.10f is true\nn = %d\n',a,b,I,err,n);
else
    fprintf('INT(sin(x),x,%f,%f) = %1.10f\nerr < %1.10f is false\nn = %d\n',a,b,I,err,n);
end