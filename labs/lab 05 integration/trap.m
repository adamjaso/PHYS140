%trap.m
%This script determines the integral of sin(x) using the trapezoidal rule

a = input('Enter the interval start: ');
b = input('Enter the interval end:   ');
err = input('Enter the error margin:   ');

xs = linspace(a,b,n);
ys = sin(xs);
I = 2;

while abs(I-(-cos(b)+cos(a))) > err
    I = 0;
    for i=2:n
        I = I + .5 * ((xs(i)-xs(i-1))*(ys(i)+ys(i-1)));
    end
    n = n + 100;
    xs = linspace(a,b,n);
    ys = sin(xs);
end

fprintf('INT(sin(x),x,%f,%f) = %1.10f\nerr = %1.10f\nn = %d\n',a,b,I,err,n);
