%% Testkatalog sammenligning af to middelværdier (kendt varians)
mu1 = 3;
n1  = 20;
mu2 = 4;
n2  = 10;
sigma = 1;
u   = 1.96;

% Data
x1 = randn(1,n1)*sigma + mu1;
x2 = randn(1,n2)*sigma + mu2;

% Parameterskøn
x1_bar = mean(x1);
x2_bar = mean(x2);

% Teststørrelse (H: mu1 = mu2)
z = (x1_bar-x2_bar)/(sigma*sqrt(1/n1+1/n2))
pval = 2*(1-normcdf(abs(z)))

% 95% konfidensinterval
d_nedre = x1_bar-x2_bar - u*sigma*sqrt(1/n1+1/n2)
d_oevre = x1_bar-x2_bar + u*sigma*sqrt(1/n1+1/n2)

figure
plot(1,x1,'b.',2,x2,'r.',1,x1_bar,'b+',2,x2_bar,'r+')
axis([0.5 2.5 1.1*min([x1 x2]) 1.1*max([x1 x2])])

%% Testkatalog sammenligning af to middelværdier (ukendt varians)
mu1 = 3;
n1  = 20;
mu2 = 4;
n2  = 10;
sigma = 1;

% Data
x1 = randn(1,n1)*sigma + mu1;
x2 = randn(1,n2)*sigma + mu2;

% Parameterskøn
x1_bar  = mean(x1);
x2_bar  = mean(x2);
s2_1    = var(x1)
s2_2    = var(x2)
s2      = 1/(n1+n2-2)*((n1-1)*s2_1+(n2-1)*s2_2)
s       = sqrt(s2)

% Teststørrelse (H: mu1 = mu2)
t = (x1_bar-x2_bar)/(s*sqrt(1/n1+1/n2))
pval = 2*(1-tcdf(abs(t),n1+n2-2))

% 95% konfidensinterval
t0 = tinv(0.975,n1+n2-2)
d_nedre = x1_bar-x2_bar - t0*s*sqrt(1/n1+1/n2)
d_oevre = x1_bar-x2_bar + t0*s*sqrt(1/n1+1/n2)

figure
plot(1,x1,'b.',2,x2,'r.',1,x1_bar,'b+',2,x2_bar,'r+')
axis([0.5 2.5 1.1*min([x1 x2]) 1.1*max([x1 x2])])
