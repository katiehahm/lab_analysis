function new_clip = wiener_clipL_adjust(largeL, currL, s, filt_pcbD, arrive_idx, last_idx, noise_thresh)
% 5/17/22
% used to adjust clip length for overlapping impacts for wiener filtering
% expands clip to data before and after, then looks at added data and 
% makes any value above noise thresh = noise thresh

diff = largeL - currL;
new_clip = filt_pcbD(arrive_idx - round(diff/2):last_idx + diff - round(diff/2),s);
make_noise_idx = find(new_clip(1:round(diff/2)) > noise_thresh);
new_clip(make_noise_idx) = noise_thresh;
make_noise_idx = find(new_clip(last_idx:end) > noise_thresh);
new_clip(make_noise_idx + last_idx - 1) = noise_thresh;

end

