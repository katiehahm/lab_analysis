% no need to edit csv file it handles headers
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Alex 4\fsr_both_regular2';
T = readtable(filepath, 'PreserveVariableNames',true);
A = table2array(T);

% columns numbers of fsr modules
person1L = 1;
person1R = 11;
trigger = 21;
person2L = 22;
person2R = 32;
personarr = [person1L, person1R, person2L, person2R];

% assign fsr variables
fsr_length = length(A(~isnan(A(:,person1L)),person1L));
fsrData = zeros(fsr_length,4);
for i = 1:4
    person = personarr(i);
    fsrData(:,i) = A(~isnan(A(:,person)),person);
end

% assign acc variables
acc_length = length(A(~isnan(A(:,(person1L + 4))),(person1L + 4)));
accData = zeros(acc_length,12);
for i = 1:4
    person = personarr(i) + 4;
    for j = 1:3
        currperson = person + j-1;
        accData(:,(i-1)*3+j) = A(~isnan(A(:,currperson)),currperson);
    end
end

trigger_signal = A(:,trigger);

Fs_fsr = 296.2963;
Fs_acc = 148.1481;
Fs_trigger = 2222.2222;

% sanity check
figure; plot(fsrData(:,1))
figure; plot(accData(:,1))
figure; plot(trigger_signal);

save(filepath,'fsrData','Fs_fsr','accData','Fs_acc','trigger_signal','Fs_trigger')