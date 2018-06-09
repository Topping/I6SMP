%% Fordelingsfunktion(CDF) til tæthedsfunktion(PDF)
clc
%Værdier
syms x u y k
%nedre grænse
a = -1; %<-----------------------------SKRIV HER
%øvre grænse
b = 1; %<------------------------------SKRIV HER
%funktion
F = 1/8*(3+4*x+x^2); %<----------------SKRIV HER
%Udregning
f = diff(F)

%% ===============Sandsynligheder===============
h = symfun(F, [x]);
syms xStor1 xMindreLigmed storXstor
xStor1 = 1-h(0) %<-------------------SKRIV HER
%xMindreLigmed = h(2/3) %<---------------SKRIV HER
%storXstor = h(0.5)-h(-0.5) %<---------------SKRIV HER

%% ===============Middelværdi===============
syms middel
middel = int(x*f,x,a,b)

%% ===============Varians===============
syms var msq
msq = int(x^2*f,x,a,b);
var = msq-middel^2

%% ===============Transformation===============
syms y
eq = y == 1/2*(x+1) %<---------------SKRIV HER
f_x = symfun(f, [x]);
% Invers
invers = solve(eq,x)
% Differentier:
differentialet = diff(invers,'y')
%grænser
nedre = solve(a == invers,y)
ovre = solve(b == invers,y)
% Formel
f_y = abs(differentialet)*f_x(invers)
expanded_f_y = expand(f_y)