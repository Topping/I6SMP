%Binomialfordeling
n = 100; % Antal trials
p = 0.5; % Sansynlighed for succes
k = 0:n; % k = antal succeser
pn = binopdf(k,n,p);

figure
bar(0:100,pn)
axis([0 100 0 0.1])
xlabel('k')
ylabel('P_n(K)')

% PDF og CDF for binomialfordeling
%   PDF: Probability Density Function (tæthedsfunktion)
%   CDF: Cumulative Distribution Function (fordelingsfunktion)
x = k;
fx = pn;
Fx = cumsum(fx);
figure
subplot(2,2,1)
plot(x,fx), title('f_X(x)')
subplot(2,2,2), title('F_X(x)')
plot(x,Fx), title('F_X(x)')