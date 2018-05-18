%trap.m
%This script determines the integral of sin(x) using the trapezoidal rule

a = input('Enter the interval start: ');
b = input('Enter the interval end:   ');
%n = input('Enter the intervals: ');

n = 1;
xs = linspace(a,b,n);
ys = sin(xs);
dx = (b-a)/n;
I = 1/2*(sin(a)+sin(b));

d1dxa=cos(a);
d1dxb=cos(b);
d3dxa=-sin(a);
d3dxb=-sin(b);

%Calculate error coefficients
c1=1/12*(d1dxa-d1dxb);
c2=1/720*(d3dxa-d3dxb);

while n <= 2^20
    %Composite Trapezoidal Method of Integration: Summation
    for i=2:n
        I = I + .5 * ((xs(i)-xs(i-1))*(ys(i)+ys(i-1)));
    end
    %Calculate error
    err=dx^2*c1-dx^4*c2;
    %Print result
    fprintf('+-------------------------------------------------------------\n');
    fprintf('|\tINT(sin(x),x,%f,%f) = %1.10f\n',a,b,I);
    fprintf('|\tn = %d\terr = %1.10f\n',n,err);
    
    %Prepare for next iteration
    n = n*2;
    xs = linspace(a,b,n);
    ys = sin(xs);
    dx = (b-a)/n;
    I = 0;
end

