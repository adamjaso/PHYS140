%simp.m
%This script determines the interval of sin(x) using Simpson's rule

a = input('Enter the interval start: ');
b = input('Enter the interval end:   ');
err = input('Enter the error margin:   ');

xs = linspace(a,b,n);

I=0;
while abs(I-(-cos(b)+cos(a))) > err
    I = 0;
    for i=1:(n-2)/2
        xis = xs(2*i:2*i+2);
        yis = sin(xis);
        x3=xis(3)^3-xis(1)^3;
        x2=xis(3)^2-xis(1)^2;
        x1=xis(3)-xis(1);
        I = I + sum([1/3 1/2 1] .* fit(xis,yis)'.*[x3 x2 x1]);
    end
    n = n + 100;
    xs = linspace(a,b,n);
    ys = sin(xs);
end

fprintf('INT(sin(x),x,%f,%f) = %1.10f\nerr = %1.10f\nn = %d\n',a,b,I,err,n);