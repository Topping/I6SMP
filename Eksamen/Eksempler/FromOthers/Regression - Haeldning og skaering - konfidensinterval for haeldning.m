%% Reggression
clc
% x og y observationer
x = [ 35 28 32 27 37 38 34 35 34 31 ]; %<---------------------------------Skriv her
y = [ 21.4 15.2 15.6 12.9 25.4 26.6 21.0 20.4 21.3 15.0 ]; %<-----Skriv her

% Udregning af hældning og skæring med matlab funktion
[r, haeldning, skaering] = regression(x, y);

% Udregnng af hældning(Beta) SIDE 435 I SMALEDE SLIDES(Model fitting)
beta = (sum((x-mean(x)).*(y-mean(y))))/(sum((x-mean(x)).^2))

% Udregning af skæring(Alpha)
alpha = mean(y)-beta*mean(x)

figure(1)
plot(x,y,'.',x,alpha+beta*x)
xlabel('X = power')
ylabel('Y = n.o. atoms')
%% Begrund at de er normalfordelte
figure(2)
qqplot(x)
figure(3)
qqplot(y)
%% Betingelser for linjær reggression opfyldt?
%Residualtegning
res = y-alpha-beta*x;
figure(4)
plot(x,res,'.',[min(x) max(x)],[0,0])
xlabel('X')
ylabel('Residuals')
%Q-Q plot
figure(5)
qqplot(res)
%% Beregn middelværdi og varians
xmean = mean(x)
xvar = var(x)
ymean = mean(y)
yvar = var(y)
%% 95% konfidensinterval for hældning
n = length(x);
% Empirical variance SIDE 436 I SAMLEDE SLIDES
sr2 = 1/(n-2)*sum((y-(skaering+haeldning*x)).^2)
sr=sqrt(sr2);
t0 = tinv(1-0.05/2,n-2)
% Beta minus og beta plus udregnes
beta_minus = haeldning - t0*sr*sqrt(1/sum((x-mean(x)).^2))
beta_plus = haeldning + t0*sr*sqrt(1/sum((x-mean(x)).^2))

%% 95% konfidens interval for Skæring - ==========IKKE TESTET==========
n = length(x)
% Empirical variance SIDE 436 I SAMLEDE SLIDES
sr2 = 1/(n-2)*sum((y-(skaering+haeldning*x)).^2);
sr=sqrt(sr2);
t0 = tinv(1-0.05/2,n-2);
alpha_minus = mean(alpha)-t0*sr*sqrt(1/n+(mean(x)^2)/(sum(x-mean(x)^2)))
alpha_minus = mean(alpha)-t0*sr*sqrt(1/n+(mean(x)^2)/(sum(x-mean(x)^2)))