function clip_energy = fft_energy(curr_window,first_clip,sec_clip,third_clip,commonF,Fs)
% 3/8/22
% takes fft of each clip, interpolates to make them same length, normalize,
% take energy and return matrix of energy

clip_energy = [];

% curr window
fft_clip = fft(curr_window);
L_clip = length(curr_window);
P2_clip = abs(fft_clip/L_clip);
P1_clip = P2_clip(1:round(L_clip/2)+1);
P1_clip(2:end-1) = 2*P1_clip(2:end-1);

f_clip = Fs*(0:(L_clip/2))/L_clip;
cutoff = findTindex(2000,f_clip);
P1_common = interp1(f_clip(1:cutoff),P1_clip(1:cutoff),commonF);
P1_norm = P1_common./max(P1_common);
clip_energy(end+1) = sum(P1_norm);

% 1st clip
if ~isempty(first_clip)
    fft_clip = fft(first_clip);
    L_clip = length(first_clip);
    P2_clip = abs(fft_clip/L_clip);
    P1_clip = P2_clip(1:round(L_clip/2)+1);
    P1_clip(2:end-1) = 2*P1_clip(2:end-1);

    f_clip = Fs*(0:(L_clip/2))/L_clip;
    cutoff = findTindex(2000,f_clip);
    P1_common = interp1(f_clip(1:cutoff),P1_clip(1:cutoff),commonF);
    P1_norm = P1_common./max(P1_common);
    clip_energy(end+1) = sum(P1_norm);
else
    clip_energy(end+1) = 0;
end

% 2nd clip
if ~isempty(sec_clip)
    fft_clip = fft(sec_clip);
    L_clip = length(sec_clip);
    P2_clip = abs(fft_clip/L_clip);
    P1_clip = P2_clip(1:round(L_clip/2)+1);
    P1_clip(2:end-1) = 2*P1_clip(2:end-1);

    f_clip = Fs*(0:(L_clip/2))/L_clip;
    cutoff = findTindex(2000,f_clip);
    P1_common = interp1(f_clip(1:cutoff),P1_clip(1:cutoff),commonF);
    P1_norm = P1_common./max(P1_common);
    clip_energy(end+1) = sum(P1_norm);
else
    clip_energy(end+1) = 0;
end

% 3rd clip
if ~isempty(third_clip)
    fft_clip = fft(third_clip);
    L_clip = length(third_clip);
    P2_clip = abs(fft_clip/L_clip);
    P1_clip = P2_clip(1:round(L_clip/2)+1);
    P1_clip(2:end-1) = 2*P1_clip(2:end-1);

    f_clip = Fs*(0:(L_clip/2))/L_clip;
    cutoff = findTindex(2000,f_clip);
    P1_common = interp1(f_clip(1:cutoff),P1_clip(1:cutoff),commonF);
    P1_norm = P1_common./max(P1_common);
    clip_energy(end+1) = sum(P1_norm);
else
    clip_energy(end+1) = 0;
end


end

