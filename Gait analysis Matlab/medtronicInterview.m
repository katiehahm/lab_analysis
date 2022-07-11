%% Question 2: ######################################
% acquire data
filename = 'data_Q2.csv';
T = readtable(filename);
A = table2array(T);

% initialize bins
numBins = 10;
minVal = min(A);
minVal = 2;
% maxVal = max(A);
binSize = (maxVal - minVal)/numBins;

% sort into bins
binCounts = zeros(1,numBins);
binMeans = zeros(1,numBins);
prevBin = minVal;
for i = 1:numBins
    binMax = prevBin + binSize;
    if i == numBins
        idx = find(A >= prevBin & A <= binMax);
    else
        idx = find(A >= prevBin & A < binMax);
    end
    binCounts(i) = length(idx);
%     binMeans(i) = (max(A(idx)) + min(A(idx)))/2;
    binMeans(i) = (prevBin + binMax)/2;
    prevBin = binMax;
end

% plot
figure;
for i = 1:numBins
    yValues = linspace(1,binCounts(i),binCounts(i));
    plot(binMeans(i),yValues,'b.')
    hold on
end
histogram(A,10)

%% Question 1: #######################################
%% read data
filename = 'data_Q1.csv';
T = readtable(filename);
A = table2array(T);

% separate train/test datasets
X = A(:,1:5);
Y = A(:,6);
N = length(Y);


N = length(Y);
Xtest = X(1:floor(N/5),:);
Xtrain = X(floor(N/5)+1:end,:);
% standardize data
XstandTrain = zeros(length(Xtrain),3);
XstandTest = zeros(floor(N/5),3);
for i = 2:4
    colData = Xtrain(:,i);
    minVal = min(colData);
    colDataSub = colData - minVal;
    maxVal = max(colDataSub);
    XstandTrain(:,i-1) = colDataSub./maxVal;
end
for i = 2:4
    colData = Xtest(:,i);
    minVal = min(colData);
    colDataSub = colData - minVal;
    maxVal = max(colDataSub);
    XstandTest(:,i-1) = colDataSub./maxVal;
end

Ytest = Y(1:floor(N/5));
Ytrain = Y(floor(N/5)+1:end);

% linear regression train
mdl = fitlm(Xtrain(:,2:end),Ytrain);
ypred = predict(mdl,Xtest(:,2:end));
% mdl = fitlm(XstandTrain,Ytrain);
% ypred = predict(mdl,XstandTest);

% evaluate
rmse = sqrt(mean((ypred - Ytest).^2));

%% plot
figure;
plot(ypred,Ytest,'r.')
xlim([5 45])
ylim([5 45])
xlabel('Predicted MEDV')
ylabel('Measured MEDV')

%% correlation between features
pcoef = corrcoef(X)
