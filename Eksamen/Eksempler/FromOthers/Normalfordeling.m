% PDF og CDF for standard normalfordeling:
%   Middelværdi = 0
%   Varians = 1
%   PDF: Probability Density Function (tæthedsfunktion)
%   CDF: Cumulative Distribution Function (fordelingsfunktion)

x = -5:0.1:5;
mean        = 0;
variance    = 1;
stddev      = sqrt(variance);   % Standard deviation
fx = normpdf(x,mean,stddev);
Fx = normcdf(x,mean,stddev);
figure
subplot(2,2,1)
plot(x,fx), title('f_X(x)')
subplot(2,2,2)
plot(x,Fx), title('F_X(x)')

% PDF og CDF for normalfordeling med
%   Middelværdi = 1
%   Varians = 2
%   PDF: Probability Density Function (tæthedsfunktion)
%   CDF: Cumulative Distribution Function (fordelingsfunktion)

x = -5:0.1:5;
mean        = 1;
variance    = 2;
stddev      = sqrt(variance);   % Standard deviation
fx = normpdf(x,mean,stddev);
Fx = normcdf(x,mean,stddev);
subplot(2,2,3)
plot(x,fx), title('f_X(x)')
subplot(2,2,4)
plot(x,Fx), title('F_X(x)')
