%% Testkatalog for uparret sammenligning
data = [
    8.0 5.6 2.4
    8.4 7.4 1.0
    8.0 7.3 0.7
    6.4 6.4 0.0
    8.6 7.5 1.1
    7.7 6.1 1.6
    7.7 6.6 1.1
    5.6 6.0 -0.4
    5.6 5.5 0.1
    6.2 5.5 0.7 ];
x1 = data(:,1);
x2 = data(:,2);

n1 = length(x1);
n2 = length(x2);

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


%% Testkatalog for parret sammenligning
data = [
    8.0 5.6 2.4
    8.4 7.4 1.0
    8.0 7.3 0.7
    6.4 6.4 0.0
    8.6 7.5 1.1
    7.7 6.1 1.6
    7.7 6.6 1.1
    5.6 6.0 -0.4
    5.6 5.5 0.1
    6.2 5.5 0.7 ];
x1 = data(:,1);
x2 = data(:,2);
d = x1-x2;

n = length(d);

% Parameterskøn
d_bar   = mean(d)
sd2     = var(d)
sd      = sqrt(s2)

% Teststørrelse (H: mu1 = mu2)
t = d_bar/(sd*sqrt(1/n))
pval = 2*(1-tcdf(abs(t),n-1))

% 95% konfidensinterval
t0 = tinv(0.975,n-1)
d_nedre = d_bar - t0*sd*sqrt(1/n)
d_oevre = d_bar + t0*sd*sqrt(1/n)