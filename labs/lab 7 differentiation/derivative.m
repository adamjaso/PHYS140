function []=derivative(x)
%derivative.m
%This function calculates the derivative of x*e^x using backward, forward
%and central differencing.

%define the h differences
d=0.005:-0.0005:0.0005;

%1ST DERIVATIVE
%calculate the derivatives for the three methods
db=dback(x,d);
df=dforw(x,d);
dc=dcent(x,d);
dc2=dcent(x,d);

%2ND DERIVATIVE
%calculate the 2nd derivatives for the three methods
d2b=d2back(x,d);
d2f=d2forw(x,d);
d2c=d2cent(x,d);

%1ST DERIVATIVE
%display the differences with their corresponding derivative approximations
fprintf('1st Derivatives\n');
for i=1:length(d)
    fprintf('h = %f   backward: df/dx = %f   forward:  df/dx = %f   central:  df/dx = %f   central2:  df/dx = %f\n',d(i),db(i),df(i),dc(i),dc2(i));
end

%2ND DERIVATIVE
%display the differences with their corresponding 2nd derivative
%approximations
fprintf('\n2nd Derivatives\n');
for i=1:length(d)
    fprintf('h = %f   backward: d2f/dx = %f   forward:  d2f/dx = %f   central:  d2f/dx = %f\n',d(i),d2b(i),d2f(i),d2c(i));
end

%1ST DERIVATIVE: EXACT 1ST DERIVATIVE
%calculate the exact value of the derivative
da=dfdx(x);

%2ND DERIVATIVE: EXACT 2ND DERIVATIVE
d2a=d2fdx(x);

%1ST DERIVATIVE: DIFFERENCES FROM EXACT
%calculate the differences of the respective approximations to the exact
%derivative
diffback=abs(db-da);
diffforw=abs(df-da);
diffcent=abs(dc-da);
%diffcent2=abs(dc2-da);

%2ND DERIVATIVE: DIFFERENCES FROM EXACT
%calculate the differences of the respective approximations to the exact
diff2back=abs(d2b-d2a);
diff2forw=abs(d2f-d2a);
diff2cent=abs(d2c-d2a);

%1ST DERIVATIVE: CHART
%plot the three differences of approximations
figure(1);
plot(d,diffback,'r+');
ylabel('| df_a_p_x/dx - df/dx |');
xlabel('h difference');
hold on;
plot(d,diffforw,'g+');
plot(d,diffcent,'b+');
title({strcat(' Error in df_a_p_x/dx at x= ',num2str(x)),'f(x)= x e^x'});

set(legend('Error in backward differencing','Error in forward differencing','Error in central differencing'),'Location','SouthEast');

%set(gca,'YScale','log');

%2ND DERIVATIVE: CHART
%plot the three differences of approximations
figure(2);
plot(d,diff2back,'r+');
ylabel('| d^2f_a_p_x/dx - df/dx |');
xlabel('h difference');
hold on;
plot(d,diff2forw,'g+');
plot(d,diff2cent,'b+');
title({strcat(' Error in d^2f_a_p_x/dx at x = ',num2str(x)),'f(x)= x e^x'});

set(legend('Error in backward differencing','Error in forward differencing','Error in central differencing'),'Location','SouthEast');


%ACCELERATION
t0=15;
t=16;
t1=20;
vt0=362.78;
vt1=517.35;

fprintf('acceleration = %f m/s^2\n',(vt1-vt0)/(t1-t0));

end

%backward difference approximation
function db=dback(x,d)
db=(f(x)-f(x-d))./d;
end
function d2b=d2back(x,d)
d2b=(f(x-2*d)-2*f(x-d)+f(x))./d.^2;
end

%forward difference approximation
function df=dforw(x,d)
df=(f(x+d)-f(x))./d;
end
function d2f=d2forw(x,d)
d2f=(f(x+2*d)-2*f(x+d)+f(x))./d.^2;
end

%central difference approximation #1
function dc=dcent(x,d)
dc=(f(x+d)-f(x-d))./(2*d);
end
function d2c=d2cent(x,d)
d2c=(f(x+d)-2*f(x)+f(x-d))./d.^2;
end

%central difference approximation #2
function dc=dcent2(x,d)
dx=d/2;
dc=(f(x+dx)-f(x-dx))./d;
end

%function definition of x*e^x
function f=f(x)
f=x.*exp(x);
end

%derivative of the function x*e^x
function dfdx=dfdx(x)
dfdx=(x+1).*exp(x);
end

%2nd derivative of the function x*e^x
function d2fdx=d2fdx(x)
d2fdx=(x+2).*exp(x);
end

%for calculation of the acceleration
function def=defderiv(x1,y1,x2,y2)
def=(y2-y1)/(x2-x1);
end