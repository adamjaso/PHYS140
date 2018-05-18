%SumOfSeries.m

terms = 10;
sinx = 0;
sum = 0;

x = 5*pi/6;

for n=1:terms
    sum = sum + (-1)^n*n/2^n;
    
    sinx = sinx + ((-1)^n*x^(2*n+1))/factorial(2*n+1);
end

fprintf('part a = %f\n',sum);
fprintf('part b = %f\n',sinx);