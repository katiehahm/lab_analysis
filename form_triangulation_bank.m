%% form delay comparison bank
% multitaps_all_A4F6 has two tap missing!!!!!
% Adjust by making first tap = 3 by copy pasting first hit
% this data is created by peaks
diff = [diff(1,:); diff(1,:); diff];
nSection = 18; % 18 because A4-F6 has 18 sections
nHits = 3; % hit each section 3 times in data
peaksBank = zeros(nSection, 5); 
for i = 0:(nSection-1)
    sums = zeros(1,5);
    for j = 1:3
        row_idx = i*nHits + j; 
        for k = 1:4 % number of sensors
            sums(k) = sums(k) + diff(row_idx, k);
        end
    end
    avgs = sums./nHits;
    avgs(5) = i+1; % number of section (1 = A4, 2 = A5, ... 18 = F6)
    peaksBank(i+1,:) = avgs;
end
save('peaksBank', 'peaksBank')