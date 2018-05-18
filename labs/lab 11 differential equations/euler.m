function []=euler(dt)
%euler.m
%This function finds the numerical solution to the ODE: y' = y^2 + 1

figure(1);
eulersimp(dt);
figure(2);
eulermod(dt);
end
