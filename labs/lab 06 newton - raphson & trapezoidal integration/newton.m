%newton.m
%Solution of x - cos(x) = 0 by the Newton-Raphson method.

%Set the initial conditions

%Correction
x = .1*pi;
dx = -(xcosx(x))/(1+sin(x));
%Derivatives
dfdx=1+sin(x);
d2fdx= cos(x);
%Error Calculation
error=-1/2*dx^2*d2fdx/dfdx; %=x-(x+dx)
%Tolerance
err = 0.000000001;
iter = 0;


%Loop through the function value
while abs(xcosx(x)) > err
    %Correction
    x = x + dx;
    dx = -(xcosx(x))/(1+sin(x));
    %Derivatives
    dfdx=1+sin(x);
    d2fdx= cos(x);
    %Error Calculation
    error=-1/2*dx^2*d2fdx/dfdx;
    fprintf('n = %d\tError = %1.40f\n',iter,error);
    %Iteration
    iter = iter+1;
end

%Return the zero of the function and the number of iterations completed
fprintf('zero = %1.9f +/- %1.40f\niterations = %d\n',x,abs(error),iter);
