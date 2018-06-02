%% Eksempel 1 - Mendels eksperiment
x = 152;
n = 580;
p0 = 1/4;
u = 1.96;

% Hypotesetest (approksimativ p-værdi)
z = (x-n*p0)/sqrt(n*p0*(1-p0))
pval = 2*(1-normcdf(abs(z)))

% Parameterskøn
p_est = x/n

% 95% konfidensinterval
p_nedre = 1/(n+u^2) * (x + u^2/2 - u*sqrt( x*(n-x)/n + u^2/4 ) )
p_oevre = 1/(n+u^2) * (x + u^2/2 + u*sqrt( x*(n-x)/n + u^2/4 ) )


%% Eksempel 2 - Drenge- og pigefødsler
x = 108;
n = 231;
p0 = 1/2;
u = 1.96;

% Hypotesetest (approksimativ p-værdi)
z = (x-n*p0)/sqrt(n*p0*(1-p0))
pval = 2*(1-normcdf(abs(z)))

% Parameterskøn
p_est = x/n

% 95% konfidensinterval
p_nedre = 1/(n+u^2) * (x + u^2/2 - u*sqrt( x*(n-x)/n + u^2/4 ) )
p_oevre = 1/(n+u^2) * (x + u^2/2 + u*sqrt( x*(n-x)/n + u^2/4 ) )

