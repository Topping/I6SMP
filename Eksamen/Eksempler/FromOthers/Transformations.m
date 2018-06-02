clc
%%
% Opgave 1 med konstant fra pdf til cdf
syms x k u
f = k*u^2;
a = -1;
b = x;
cdf = int(f,u,a,b)
clear cdf
% Fra cdf til pdf

cdf = 1/8*(3+4*x+x^2);
pdf = diff(cdf)

%%
% Transformationssætning 
clear
% f_y(y) er tæthedsfunktionen for den stokastiske variabel Y givet ved Y=
% arcsin(X)

% Vis at f_y(y) er givet ved : something something...

syms X y x 
eq = y == 1/4*(x+1/2);

f_x(x) = 3/2 * x^2;
% X interval grænser a <= x <= b
a = -1;
b = 1;

% Invers
invers = solve(eq,x)

% Differentier:

differentialet = diff(invers,'y')


 nedre = solve(a == invers,y)
 ovre = solve(b == invers,y)
% Formel


f_y = abs(differentialet)*f_x(invers)
pdf_res = expand(f_y)


