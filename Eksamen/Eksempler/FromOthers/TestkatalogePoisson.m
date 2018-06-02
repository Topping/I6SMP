%% Eksempel 3 - Forskellige Poisson fordelinger (t�thedsfunktioner)
x = 1:30;
fx1 = poisspdf(x,3);
fx2 = poisspdf(x,6);
fx3 = poisspdf(x,9);
fx4 = poisspdf(x,12);
fx5 = poisspdf(x,15);
figure
plot(x,fx1,x,fx2,x,fx3,x,fx4,x,fx5)
xlabel('x = antal ankomster per tid')
ylabel('Sandsynlighed')
legend('\lambda=3','\lambda=6','\lambda=9','\lambda=12','\lambda=15')

%% Eksempel 4 - Anvendelse af Poisson fordelingen
% En butik har 300 bes�gende p� 2 timer. Hvad er sandsynligheden for, at
% der kommer mere end 170 bes�gende den n�ste time? 
% Data:             x=300
% Parametersk�n:	lambda = 300/2 = 150 bes�gende/time
xobs = 300
t = 2
lambda_est = xobs/t
x = 170
prob = 1-poisscdf(x,lambda_est)

%% Eksempel 5 - Radioaktive henfald
data = [
    0 57
    1 203
    2 383
    3 525
    4 532
    5 408
    6 273
    7 139
    8 45
    9 27
    10 10
    11 4
    12 0
    13 1
    14 1 ];

t = 72;
t_total = t*sum(data(:,2))    % = 72*2806
x = sum(data(1,:).*data(2,:))
u = 1.96; % 95% konfidenceinterval
% Parametersk�n
lambda_est = x/t_total

% 95% konfidensinterval
lambda_nedre = 1/t_total * (x + u^2/2 - u*sqrt( x + u^2/4 ) )
lambda_oevre = 1/t_total * (x + u^2/2 + u*sqrt( x + u^2/4 ) )