%newton.m
%Solution of x - cos(x) = 0 by the Newton-Raphson method.

%Set the initial conditions
x = .1*pi;
dx = -(xcosx(x))/(1+sin(x));
err = 0.000000001;
iter = 0;

%Loop through the function value
while abs(xcosx(x)) > err
    x = x + dx;
    dx = -(xcosx(x))/(1+sin(x));
    iter = iter+1;
end

%Return the zero of the function and the number of iterations completed
fprintf('zero = %1.9f\niterations = %d\n',x,iter);
