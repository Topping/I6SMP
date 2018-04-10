%% ETSMP Ebola Outbreak II
% Udvalgte løsninger

%% Process 1:
% y(n)=x

%% Creations of realisations
days=10;
patients=10;
propability_success=0.5;
yn=binornd(patients,propability_success,1); %binomial random variables

%% Ensemble mean, brugt side 48 i formelsamling
% For binomial: E[x]=n*p
Ensemble_mean=patients*propability_success;  %= 5, s 48 i formelsamling

%% Ensemble variance, brugt side 48 i formelsamling
% var(x) = p*n(1-p)
Ensemble_variance=patients*propability_success*(1-propability_success); %=2.5

%% Verifikation med matlab
yn=binornd(patients,propability_success,1,100000);
mean(yn); %Should be equal to the ensemble mean
var(yn); % Should be equal to the ensemble variance

%% WSS or ergodic?
%% The process is WSS as the mean and variance is constant with time.
%% The process is not ergodic, as one realization has a variance of 0.

%% Process 2
% yn=x+wn
x=binornd(patients,propability_success,1);
wn=randi([-2 2],1,days); %creates discrete uniformly distributed data.
yn=x+wn;

%% The ensemble mean and varaince?
% da E[wn]= (a+b)/2=(2+-2)/2=0, E[x]= n*p = 5
% E[yn]=E[x+wn]=E[x]+E[wn]= 0+5= 5
%
% da var(wn)=((b-a+1)^2-1)/12= ((2--2+1)^2-1)/12=2 
% var(yn)=E[yn^2]-E[yn]^2=var(x)+var(wn)=p*n(1-p)+2=4.5

%% Verifikation med matlab funktion, vi bør gøre dette for alle 10 dage, men da processen er WSS kan vi nøjes med dag 1
testDag1=binornd(patients,propability_success,1,10000)+randi([-2 2],1,10000);
mean(testDag1) %skal gerne give 5
var(testDag1) % skal gerne give 4.5

%% WSS or ergodic?
%% The process is WSS as the mean and variance is constant with time.
%% The process is not ergodic, as the timely mean is not always equal to the ensemle mean.
