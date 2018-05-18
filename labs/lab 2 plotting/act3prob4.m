T=[5840 22400 13260 9400 3130 4280 28200]/5480;
LLdata=[1 13400 150 108 0.004 0.15 3400];
RRdata=[1 7.8 3.5 3.7 0.18 0.76 8];
LLsun=(RRdata.^2).*(T.^4);
p1=plot(LLdata,T,'c*');
xlabel('Temperature Ratio (T / T_s_u_n)');
ylabel('Ratio relative to Sun (L / L_s_u_n and R / R_s_u_n)');
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'XDir','reverse');
hold on;
p2=plot(LLsun,T,'mo');
set(gca,'YScale','log');
set(gca,'XDir','reverse')
legend('L / L_s_u_n','R / R_s_u_n');
title('Hertzprung - Russell Diagram');