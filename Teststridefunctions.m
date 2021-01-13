%% slight limp
filenames = {'data/data_10-31-2020_17-08-30','data/data_10-31-2020_17-09-11','data/data_10-31-2020_17-10-07','data/data_10-31-2020_17-10-43','data/data_10-31-2020_17-11-18','data/data_10-31-2020_17-12-30',...
    'data/data_10-31-2020_17-13-45','data/data_10-31-2020_17-14-14','data/data_10-31-2020_17-14-53','data/data_10-31-2020_17-15-19','data/data_10-31-2020_17-15-55',...
    'data/data_10-31-2020_17-16-24','data/data_10-31-2020_17-17-12','data/data_10-31-2020_17-17-43','data/data_10-31-2020_17-18-11','data/data_10-31-2020_17-18-48','data/data_10-31-2020_17-19-11',...
    'data/data_10-31-2020_17-19-36','data/data_10-31-2020_17-20-04'};
impactN=6;
[diff_foot1, diff_foot2] = hist_stridesamefoot(filenames,impactN)
% [avgpeaks1,avgpeaks2,avgpeaktable]=peakbyfoot(filenames,impactN)
% [differences1, differences2, differences3, mu1, mu2,mu3,sigma1,sigma2,sigma3] = hist_stridetime(filenames,impactN)
% [peaks1, peaks2, peaks3,mu1,sigma1, mu2, sigma2,mu3,sigma3] = hist_stridemag(filenames,impactN)
%% no limp
filenames = {'data/data_10-19-2020_15-49.mat', 'data/data_10-19-2020_15-56.mat','data/data_10-19-2020_15-57.mat','data/data_10-19-2020_15-58.mat','data/data_10-19-2020_15-59.mat','data/data_10-19-2020_16-00.mat',...
    'data/data_10-19-2020_16-01.mat','data/data_10-19-2020_16-02.mat','data/data_10-19-2020_16-03.mat','data/data_10-19-2020_16-04.mat','data/data_10-19-2020_16-05.mat',...
    'data/data_10-19-2020_16-06.mat','data/data_10-19-2020_16-07.mat','data/data_10-19-2020_16-08.mat','data/data_10-19-2020_16-09.mat','data/data_10-19-2020_16-10.mat','data/data_10-19-2020_16-11.mat',...
    'data/data_10-19-2020_16-12.mat','data/data_10-19-2020_16-13.mat','data/data_10-19-2020_16-14.mat','data/data_10-19-2020_16-15.mat','data/data_10-19-2020_16-16.mat','data/data_10-19-2020_16-17.mat','data/data_10-19-2020_16-18.mat',...
    'data/data_10-19-2020_16-19.mat','data/data_10-19-2020_16-20.mat','data/data_10-19-2020_16-21.mat','data/data_10-19-2020_16-22.mat','data/data_10-19-2020_16-23.mat','data/data_10-19-2020_16-24.mat',...
    'data/data_10-19-2020_16-25.mat','data/data_10-19-2020_16-26.mat','data/data_10-19-2020_16-27.mat','data/data_10-19-2020_16-28.mat','data/data_10-19-2020_16-29.mat','data/data_10-19-2020_16-30.mat',...
    'data/data_10-19-2020_16-31.mat','data/data_10-19-2020_16-32.mat','data/data_10-19-2020_16-33.mat','data/data_10-19-2020_16-34.mat','data/data_10-19-2020_16-35.mat','data/data_10-19-2020_16-36.mat',...
    'data/data_10-19-2020_16-37.mat','data/data_10-19-2020_16-38.mat','data/data_10-19-2020_16-39.mat','data/data_10-19-2020_16-40.mat','data/data_10-19-2020_16-41.mat','data/data_10-19-2020_16-42.mat','data/data_10-19-2020_16-43.mat','data/data_10-19-2020_16-44.mat'}; 
impactN=6;
[diff_foot1, diff_foot2] = hist_stridesamefoot(filenames,impactN)
% [avgpeaks1,avgpeaks2,avgpeaktable]=peakbyfoot(filenames,impactN)
% [differences1, differences2, differences3, mu1, mu2,mu3,sigma1,sigma2,sigma3] = hist_stridetime(filenames,impactN)
%[peaks1, peaks2, peaks3,mu1,sigma1, mu2, sigma2,mu3,sigma3] = hist_stridemag(filenames,impactN)
%% severe limp
filenames ={'data/data_10-31-2020_17-22-35','data/data_10-31-2020_17-23-01','data/data_10-31-2020_17-23-26','data/data_10-31-2020_17-23-48','data/data_10-31-2020_17-24-16','data/data_10-31-2020_17-24-38',...
    'data/data_10-31-2020_17-25-04','data/data_10-31-2020_17-25-29','data/data_10-31-2020_17-25-53','data/data_10-31-2020_17-26-23','data/data_10-31-2020_17-26-48','data/data_10-31-2020_17-27-16',...
    'data/data_10-31-2020_17-27-45','data/data_10-31-2020_17-28-13','data/data_10-31-2020_17-28-45','data/data_10-31-2020_17-29-19','data/data_10-31-2020_17-29-42','data/data_10-31-2020_17-30-08',...
    'data/data_10-31-2020_17-30-36','data/data_10-31-2020_17-31-10'};
impactN=6;
[diff_foot1, diff_foot2] = hist_stridesamefoot(filenames,impactN)
% [avgpeaks1,avgpeaks2, avgpeaktable]=peakbyfoot(filenames,impactN)
% [differences1, differences2, differences3, mu1, mu2,mu3,sigma1,sigma2,sigma3] = hist_stridetime(filenames,impactN)
%[peaks1, peaks2, peaks3,mu1,sigma1, mu2, sigma2,mu3,sigma3] = hist_stridemag(filenames,impactN)
%% same speed limp?
filenames={'data/data_11-06-2020_11-58-08','data/data_11-06-2020_11-58-32','data/data_11-06-2020_11-58-54','data/data_11-06-2020_11-59-18','data/data_11-06-2020_11-59-51'}
impactN=6;
[avgpeaks1,avgpeaks2]=peakbyfoot(filenames,impactN)
%[peaks1, peaks2, peaks3,mu1,sigma1, mu2, sigma2,mu3,sigma3] = hist_stridemag(filenames,impactN)
%[differences1, differences2, differences3, mu1, mu2,mu3,sigma1,sigma2,sigma3] = hist_stridetime(filenames,impactN)
%% test distance travelled function
nolimpfilename ={'data/data_10-19-2020_15-49.mat'};%no limp
slightlimpfilename ={'data/data_10-31-2020_17-08-30.mat'}; %slight limp
severelimpfilename ={'data/data_10-31-2020_17-22-35.mat'}; %severe limp
totallength = 6; %total length walked

[strides]= distancetravelled(nolimpfilename,slightlimpfilename, severelimpfilename, totallength)
%% test stride speeds function
nolimpfilename ={'data/data_10-19-2020_15-49.mat'};%no limp
slightlimpfilename ={'data/data_10-31-2020_17-08-30'}; %slight limp
severelimpfilename ={'data/data_10-31-2020_17-22-35'}; %severe limp
impactN=6

[foot1, foot2,peak_diffandstepinfo] = stridetime_anaylsisfunction(peak_idx,impactN,trainedClassifier)
% [totalspeeds] = stridespeeds(nolimpfilename,slightlimpfilename, severelimpfilename,steplength)
% steplength = 1;%consistent step length