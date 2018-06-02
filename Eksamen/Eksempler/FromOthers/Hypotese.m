
%% Hypotese for poisson fordeling

lambda = 400;       % <---------------------- skriv lige d�r!
x = 1140;           % <---------------------- skriv lige d�r!
t = 3;              % <---------------------- skriv lige d�r!

% Test st�rrelse for poisson fordeling

z_poiss = (x-t*lambda)/sqrt(t*lambda)

%% Test st�rrelse for binomial fordeling 
n = 0; 
p0 = 0;
z_bino = (x-n*p0)/sqrt(n*p0*(1-p0))

%% Hypotese test for normalfordeling med kendt varians

x_hat   = 300          % observerede resultat
my      = 310          % p�st�ede / 0-hypotese
sigma   = 50           % standard afvigelse
n       = 40           % antal observationer

% Test st�rrelse
z_norm =  (x_hat - my) / (sigma/sqrt(n))



%% Udregner approx p-value, HUSK at inds�tte z-v�rdi lige herunder!
z = z_poiss;         % <---------------------- skriv lige d�r!
p_value = 2 * abs(1-normcdf(abs(z)))


%% Hypotese test ved sammenligning af to datas�t med ukendt f�lles varians

x = 0;      % datas�t 1  
y = 0;      % datas�t 2

xmean = mean(X);
xvar = var(X);
ymean = mean(Y);
yvar = var(Y);
nx = length(X);
ny = length(Y);
s = 1/(nx+ny-2)*((nx-1)*xvar + (ny-1)*yvar)  = (xmean-ymean)/(s*sqrt(1/nx+1/ny))
pval = 2*(1-tcdf(abs(t),nx+ny-2))



