%% 9/23/20 Logistic Regression Classifier
set(0,'DefaultFigureVisible','off')
drop = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-18-2020_11-22.mat'); % dropping soccer ball data
labels = zeros(10, 1);
filt_data1 = lpf_data(drop.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
impact_magnitude = sum(peak_val1, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
table1 = table(impact_magnitude, peak_width, labels);

slam = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-18-2020_11-23.mat'); % slamming soccer ball data
labels = ones(10, 1);
filt_data2 = lpf_data(slam.datas);
clean_data2 = clean_envelope(filt_data2,Fs);
[onset_idx2, peak_idx2, peak_val2] = TDOA(clean_data2,impactN,Fs,loc_names);
impact_magnitude = sum(peak_val2, 2);
peak_width = 2*abs(mean(onset_idx2,2)- mean(peak_idx2,2));
table2 = table(impact_magnitude, peak_width, labels);

final_table = vertcat(table1, table2);
disp(final_table)

%% 9/30/20 Classifier with new data @ different spots
soccer_drop_1 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-27.mat'); % BQ & CQ
labels = zeros(10, 1);
filt_data1 = lpf_data(soccer_drop_1.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(3/14);
a5 = peak_val1(:,[2])*(5/14);
e1 = peak_val1(:,[3])*(6/14);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_drop_table1 = table(impact_magnitude,peak_width, labels);

soccer_drop_2 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-28.mat'); % DQ & EQ
labels = zeros(10, 1);
filt_data1 = lpf_data(soccer_drop_2.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(6/17);
a5 = peak_val1(:,[2])*(7/17);
e1 = peak_val1(:,[3])*(4/17);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_drop_table2 = table(impact_magnitude, peak_width, labels);

soccer_drop_3 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-29.mat'); % ER & ES
labels = zeros(10, 1);
filt_data1 = lpf_data(soccer_drop_3.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(7/18);
a5 = peak_val1(:,[2])*(6/18);
e1 = peak_val1(:,[3])*(5/18);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_drop_table3 = table(impact_magnitude, peak_width, labels);

soccer_drop_table = vertcat(soccer_drop_table1, soccer_drop_table2, soccer_drop_table3);



soccer_bounce_1 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-30.mat'); % BQ & CQ
labels = ones(10, 1);
filt_data1 = lpf_data(soccer_bounce_1.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(3/14);
a5 = peak_val1(:,[2])*(5/14);
e1 = peak_val1(:,[3])*(6/14);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_bounce_table1 = table(impact_magnitude, peak_width, labels);

soccer_bounce_2 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-33.mat'); % DQ & EQ
labels = ones(10, 1);
filt_data1 = lpf_data(soccer_bounce_2.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(6/17);
a5 = peak_val1(:,[2])*(7/17);
e1 = peak_val1(:,[3])*(4/17);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_bounce_table2 = table(impact_magnitude, peak_width, labels);

soccer_bounce_3 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-31.mat'); % ER & ES
labels = ones(10, 1);
filt_data1 = lpf_data(soccer_bounce_3.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(7/18);
a5 = peak_val1(:,[2])*(6/18);
e1 = peak_val1(:,[3])*(5/18);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_bounce_table3 = table(impact_magnitude, peak_width, labels);

soccer_bounce_table = vertcat(soccer_bounce_table1, soccer_bounce_table2, soccer_bounce_table3);

soccer_table = vertcat(soccer_drop_table, soccer_bounce_table)

%% 9/30/2020 Box vs Ball Drop
soccer_drop_1 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-27.mat'); % BQ & CQ
labels = zeros(10, 1);
filt_data1 = lpf_data(soccer_drop_1.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(3/14);
a5 = peak_val1(:,[2])*(5/14);
e1 = peak_val1(:,[3])*(6/14);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_drop_table1 = table(impact_magnitude,peak_width, labels);

soccer_drop_2 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-28.mat'); % DQ & EQ
labels = zeros(10, 1);
filt_data1 = lpf_data(soccer_drop_2.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(6/17);
a5 = peak_val1(:,[2])*(7/17);
e1 = peak_val1(:,[3])*(4/17);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_drop_table2 = table(impact_magnitude, peak_width, labels);

soccer_drop_3 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-29.mat'); % ER & ES
labels = zeros(10, 1);
filt_data1 = lpf_data(soccer_drop_3.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(7/18);
a5 = peak_val1(:,[2])*(6/18);
e1 = peak_val1(:,[3])*(5/18);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
soccer_drop_table3 = table(impact_magnitude, peak_width, labels);

soccer_drop_table = vertcat(soccer_drop_table1, soccer_drop_table2, soccer_drop_table3);


box_drop_1 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-38.mat'); % BQ & CQ
labels = ones(10, 1);
filt_data1 = lpf_data(box_drop_1.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(3/14);
a5 = peak_val1(:,[2])*(5/14);
e1 = peak_val1(:,[3])*(6/14);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
bobx_drop_table1 = table(impact_magnitude,peak_width, labels);


box_drop_3 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-52.mat'); % ER & ES
labels = ones(10, 1);
filt_data1 = lpf_data(box_drop_3.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 10; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(7/18);
a5 = peak_val1(:,[2])*(6/18);
e1 = peak_val1(:,[3])*(5/18);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
box_drop_table3 = table(impact_magnitude, peak_width, labels);

box_drop_table = vertcat(box_drop_table1, box_drop_table3);

final_drop = vertcat(soccer_drop_table, box_drop_table)


%% 9/30/2020
steps_1 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-55.mat'); % BQ & CQ
labels = zeros(6, 1);
filt_data1 = lpf_data(steps_1.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 6; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(3/14);
a5 = peak_val1(:,[2])*(5/14);
e1 = peak_val1(:,[3])*(6/14);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
steps_table1 = table(impact_magnitude,peak_width, labels);

steps_2 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-56.mat'); % DQ & EQ
labels = zeros(6, 1);
filt_data1 = lpf_data(steps_2.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 6; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(6/17);
a5 = peak_val1(:,[2])*(7/17);
e1 = peak_val1(:,[3])*(4/17);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
steps_table2 = table(impact_magnitude, peak_width, labels);


stomps_1 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-57.mat'); % BQ & CQ
labels = ones(6, 1);
filt_data1 = lpf_data(stomps_1.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 6; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(3/14);
a5 = peak_val1(:,[2])*(5/14);
e1 = peak_val1(:,[3])*(6/14);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
stomps_table1 = table(impact_magnitude,peak_width, labels);

stomps_2 = load('C:\Users\bpdwy\Documents\MIT\Academic\2020-2021\UROP\Analysis\data\data_09-25-2020_15-58.mat'); % DQ & EQ
labels = ones(6, 1);
filt_data1 = lpf_data(stomps_2.datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
clean_data1 = clean_envelope(filt_data1,Fs);

impactN = 6; % variable input value
[onset_idx1, peak_idx1, peak_val1] = TDOA(clean_data1,impactN,Fs,loc_names);
a1 = peak_val1(:,[1])*(6/17);
a5 = peak_val1(:,[2])*(7/17);
e1 = peak_val1(:,[3])*(4/17);
peak_vals = horzcat(a1, a5, e1);
impact_magnitude = sum(peak_vals, 2);
peak_width = 2*abs(mean(onset_idx1,2)- mean(peak_idx1,2));
stomps_table2 = table(impact_magnitude, peak_width, labels);

steps_stomps = vertcat(steps_table1, steps_table2, stomps_table1, stomps_table2)

%% 10/7/20 Classifiying jumps vs stomps & steps

