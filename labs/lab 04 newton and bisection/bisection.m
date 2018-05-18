%bisection.m
%Solution of the equation x - cos(x) = 0 by the bisection method.

%Set the initial conditions
lo = .01*pi;
hi = pi/2;
m = (hi+lo)/2;
err = 0.000000001;
iter = 0;

%Recursively evaluate the function until it returns a value less than the
%error margin
while abs(xcosx(hi)) > err && abs(xcosx(m)) > err && abs(xcosx(lo)) > err
    if xcosx(hi)*xcosx(m) < 0
        lo = m;
        m = (hi+lo)/2;
    else
        hi = m;
        m = (hi+lo)/2;
    end
    iter = iter+1;
end

fprintf('zero = %1.9f\niterations = %d\n',m,iter);