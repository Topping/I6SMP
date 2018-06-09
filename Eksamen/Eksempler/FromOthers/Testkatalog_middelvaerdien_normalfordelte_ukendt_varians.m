
% Testkatalog for middelværdien i normalfordelte data (ukendt varians))

x = [ 5.36 5.29 5.58 5.65 5.57 5.53 5.62 5.29 ...
      5.44 5.34 5.79 5.10 5.27 5.39 5.42 5.47 ...
      5.63 5.34 5.46 5.30 5.75 5.68 5.85 ];

n = length(x); % antal samples
mu0 = 5.517; 
s2 = var(x) %s2 = impirisk varians

% Teststørrelse
x_hat = mean(x)
t = (x_hat-mu0)/sqrt(s2/n) 
pval = 2*(1-tcdf(abs(t),n-1)) % p-værdi

% 95% konfidensinterval
t0 = tinv(0.975,n-1)
mu_nedre = x_hat - t0*sqrt(s2/n)
mu_oevre = x_hat + t0*sqrt(s2/n)

qqplot(x)