% PDF og CDF for uniformfordeling
%   PDF: Probability Density Function (tæthedsfunktion)
%   CDF: Cumulative Distribution Function (fordelingsfunktion)

x = 0:0.1:10;
fx = zeros(length(x));
a = 2.;
b = 6.;
fx(find(x==a):find(x==b))=1;
fx = fx/sum(fx);    % Normaliser areal under kurven til 1
Fx = cumsum(fx);    % Akkumuleret sum af elementerne i fx
subplot(2,2,3)
plot(x,fx), title('f_X(x)')
subplot(2,2,4)
plot(x,Fx), title('F_X(x)')
