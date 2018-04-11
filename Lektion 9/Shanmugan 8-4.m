%Solution for Shanmugan 8.4

C1=[3.42, 3.48, 3.48, 3.54, 3.51, 3.48, 3.57, 3.59, 3.63, 3.50, 3.45, 3.51, 3.55, 3.59, 3.50, 3.61];
meanC1=sum(C1)/length(C1); %Calculation of mean
mean(C1) % verifikation

varC1=1/(length(C1)-1).*sum((C1-meanC1).^2); %Calculation of variance
var(C1)  % verifikation