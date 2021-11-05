data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\';
name = [data_root_katie,'Subject 3\Run_number_277_Plot_and_Store_Rep_3.4.csv'];
load([data_root_katie,'Subject 3\subj3_fsr_insole'])
filename = [data_root_katie,'Subject 3\subj3_fsr_insole'];
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
Data(1:end,9) = A(1:end,18);
Data = transpose(Data);

% save(filename, 'Data','Time','Channels','Fs')