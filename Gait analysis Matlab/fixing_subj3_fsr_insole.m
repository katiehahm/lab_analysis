% this file is used to correct it when fsr process didn't work to .mat
% so the file was save to .xlsx instead

filename = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\Subj 2\subj2_';
name = [filename,'fsr_regular2_csv.csv'];
load([filename,'fsr_regular2_mat'])

T = readtable(name);
A = table2array(T);
Data = zeros(length(A),9);
Data(:,1) = A(:,2);
Data(:,2) = A(:,4);
Data(:,3) = A(:,6);
Data(:,4) = A(:,8);
Data(:,5) = A(:,10);
Data(:,6) = A(:,12);
Data(:,7) = A(:,14);
Data(:,8) = A(:,16);
Data(:,9) = A(:,18);
Data = transpose(Data);

filename = [filename,'fsr_regular2'];
save(filename, 'Data','Time','Channels','Fs')