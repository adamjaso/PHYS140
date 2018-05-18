function []=interpol()
%interpol.m
%This function interpolates the function f(x) = 1.5^x cos(2x);
Xs=0:1:5;
Ys=[1.0 -.6242 -1.4707 3.2406 -.7366 -6.3717];
xs=0:.01:5;
ys=exact(xs);

N=length(xs);
ps(N)=0;
for i=1:N
    ps(i)=p(xs(i),Xs,Ys);
end

plot(Xs,Ys,'b+');
xlabel('x values');
ylabel('f(x)');
title({'Plot of the function  f(x) = 1.5^x cos(2x)','with interpolated points'});

hold on;
plot(xs,ps,'r-');
plot(xs,ys,'g-');
set(legend('Interpolated Points','Interpolated Function: p(x)','Actual Function: f(x)'),'Location','SouthWest');

fprintf('x       \tp(x)    \ty(x)    \t|p-y|   \n--------    --------    --------    --------\n');
for i=1:N
   fprintf('%f\t%f\t%f\t%f\n',xs(i),ps(i),ys(i),abs(ps(i)-ys(i)));
end

end

function p=p(x,xs,ys)
N=length(xs);
p=0;
for j=1:N
    p=p+L(x,xs,j)*ys(j);
end
end

function L=L(x,xs,j)
N=length(xs);
%fprintf('%d\n',N);
L=1;
for i=1:N
    if i ~= j
        L=L*(x-xs(i))/(xs(j)-xs(i));
    end
end
end

function vs=exact(xs)
vs=1.5.^xs.*cos(2*xs);
end