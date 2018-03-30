clc
clear
%% SAMPLINGSOPGAVE 07-03-2018 SMP

% Create random sample values
x = rand(1000,1);

% Random sample values from a Rayleigh distribution
Y1 = sampling_distribution(x, 'Rayleigh');

% Random sample values from a exponential distribution
Y2 = sampling_distribution(x, 'exp');


%% SUMOPGAVE

