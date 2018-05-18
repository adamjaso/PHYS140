%ls.m
%This program requests a set of coordinate pairs with x values in row one
%and y values in row two.

x=input('Enter your x values = ');
y=input('Enter your y values = ');

X=sum(x);
Y=sum(y);
X2=sum(x.^2);
XY=sum(x.*y);
n=length(x);

coeff=[X2,X;X,n]^(-1)*[XY;Y];

fprintf('a = %f\nb = %f\n',coeff(1),coeff(2));

plot(x,y,'ob');
hold on;
plot(x,coeff(1)*x+coeff(2),'-r');
xlabel('time');
ylabel('velocity');
title('Velocity vs. Time Data & Fitted Line');
legend('Data','Fitted Line');

eq='y = a x + b';
first=['a = ' num2str(coeff(1))];
second=['b = ' num2str(coeff(2))];
annotation('textbox',[.2 .7 .3 .15],'String',strvcat(eq,first,second));