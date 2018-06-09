%% Sandsynigheder i distributioner
%% Poisson fordeling
clc
% Tidsperiode
%t = 100; %<--------------------Skriv her
%Antal observationer
%x = 25; %<---------------------Skriv her

lambda = 3

%% Sandsynlighed
% antal i sandsynlighed X>x1
lambda = 3
x = 5;
p = poisscdf(x,lambda);
sandsynlighed = 1-p

%% antal i sandsynlighed X<=x2
x = 7; %<----------------------Skriv her
p = poisscdf(x,lambda)

%% antal i sandsynlighed X=x
x = 0; %<----------------------Skriv her
p = poisspdf(x,lambda)