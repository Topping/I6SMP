%%
clc
clear
%Tæthedsfunktion(PDF) til fordelingsfunktion(CDF)
%værdier
syms x u y k

a = -1;
b = y;
f = 3/8*(u^2+1);
F = int(f,u,a,b)

%% ===============Sandsynligheder===============
h = symfun(F, [y]); 
syms YStor1 YMindreLigmed storYstor
YStor1 = 1-h(1/4)
YMindreLigmed = h(-1/2)
%storXstor = h(0.5)-h(-0.5)

%% ===============Middelværdi===============
clear b x f
syms x middel l f

a = -1;
b = 1;
f = 2/3*y;
middel = int(x*f,y,a,b)

%% ===============Varians===============
syms var msq
msq = int(x^2*f,y,a,b)

var = msq-middel^2

%% ===============Transformation===============
syms y x 
eq = y == x^2; %<------------------SKRIV HER
f_x = symfun(f, [x]);
% Invers
invers = solve(eq,x)
% Differentier:
differentialet = diff(invers,'y')
nedre = solve(a == invers,y)
ovre = solve(b == invers,y)
% Formel
f_y = abs(differentialet)*f_x(invers)
expanded_f_y = expand(f_y)