function [d] = noise_to_zero(d,Fs)
% uses the first 1 sec of data for noise
% sets anything in this range = 0

for i = 1:4
    noise = d(1:Fs,i);
    maxv = max(noise);
    minv = min(noise);
    noisei = find(d(:,i) < maxv & d(:,i) > minv);
    d(noisei,i) = 0;
end
end

