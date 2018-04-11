%% Solution for Shanmugan problem 8.2

%Generation of data
meanX=2;
varX=9;
samplesize=100;
data=sqrt(varX).*randn(samplesize,1)+meanX;
hist(data)