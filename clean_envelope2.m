function [data] = clean_envelope2(d, Fs)
% takes the first second of data to be noise
% takes any signal within that threshold to be 0
% takes envelope of remaining data
    data = d;

    noise1 = d(1:floor(Fs),1);
    noise2 = d(1:floor(Fs),2);
    noise3 = d(1:floor(Fs),3);
    noise4 = d(1:floor(Fs),4);
    max1 = max(noise1);
    min1 = min(noise1);
    max2 = max(noise2);
    min2 = min(noise2);
    max3 = max(noise3);
    min3 = min(noise3);
    max4 = max(noise4);
    min4 = min(noise4);
    
    i1 = find(d(:,1) < max1 & d(:,1) > min1);
    d1 = d(:,1);
    d1(i1) = 0;
    data(:,1) = abs(d1);
    
    i2 = find(d(:,2) < max2 & d(:,2) > min2);
    d2 = d(:,2);
    d2(i2) = 0;
    data(:,2) = abs(d2);
    
    i3 = find(d(:,3) < max3 & d(:,3) > min3);
    d3 = d(:,3);
    d3(i3) = 0;
    data(:,3) = abs(d3);
    
     i4 = find(d(:,4) < max4 & d(:,4) > min4);
    d4 = d(:,4);
    d3(i4) = 0;
    data(:,4) = abs(d4);
    
    
    
end