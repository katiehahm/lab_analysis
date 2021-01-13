%Normalize distance
up = 11/4
across = 13/5
orig=[0,0]
distance3=[norm([0 4*up]-orig), norm([1*across 4*up]-orig), norm([2*across 4*up]-orig), norm([3*across 4*up]-orig), norm([4*across 4*up]-orig), norm([5*across 4*up]-orig);...
           norm([0 3*up]-orig), norm([1*across 3*up]-orig), norm([2*across 3*up]-orig), norm([3*across 3*up]-orig), norm([4*across 3*up]-orig), norm([5*across 3*up]-orig);...
           norm([0 2*up]-orig), norm([1*across 2*up]-orig), norm([2*across 2*up]-orig), norm([3*across 2*up]-orig), norm([4*across 2*up]-orig), norm([5*across 2*up]-orig);...
           norm([0 1*up]-orig), norm([1*across 1*up]-orig), norm([2*across 1*up]-orig), norm([3*across 1*up]-orig), norm([4*across 1*up]-orig), norm([5*across 1*up]-orig);...
           0, norm([1*across 0]-orig), norm([2*across 0]-orig), norm([3*across 0]-orig), norm([4*across 0]-orig), norm([5*across 0]-orig)];
normdistance3 = distance3 - min(distance3(:))
normdistance3 = normdistance3 ./ max(normdistance3(:)) %normalized distance from 3rd sensor
       
distance2=[norm([0 0]-orig), norm([1*across 0]-orig), norm([2*across 0]-orig), norm([3*across 0]-orig), norm([4*across 0]-orig), norm([5*across 0]-orig);...
           norm([0 1*up]-orig), norm([1*across 1*up]-orig), norm([2*across 1*up]-orig), norm([3*across 1*up]-orig), norm([4*across 1*up]-orig), norm([5*across 1*up]-orig);...
           norm([0 2*up]-orig), norm([1*across 2*up]-orig), norm([2*across 2*up]-orig), norm([3*across 2*up]-orig), norm([4*across 2*up]-orig), norm([5*across 2*up]-orig);...
           norm([0 3*up]-orig), norm([1*across 3*up]-orig), norm([2*across 3*up]-orig), norm([3*across 3*up]-orig), norm([4*across 3*up]-orig), norm([5*across 3*up]-orig);...
           norm([0 4*up]-orig), norm([1*across 4*up]-orig), norm([2*across 4*up]-orig), norm([3*across 4*up]-orig), norm([4*across 4*up]-orig), norm([5*across 4*up]-orig)];
normdistance2 = distance2 - min(distance2(:))
normdistance2 = normdistance2 ./ max(normdistance2(:)) %normalized distance from 2rd sensor
       
distance1=[norm([5*across 0]-orig), norm([4*across 0]-orig), norm([3*across 0]-orig), norm([2*across 0]-orig), norm([1*across 0]-orig), norm([0 0]-orig);...
           norm([5*across 1*up]-orig), norm([4*across 1*up]-orig), norm([3*across 1*up]-orig), norm([2*across 1*up]-orig), norm([1*across 1*up]-orig), norm([0 1*up]-orig);...
           norm([5*across 2*up]-orig), norm([4*across 2*up]-orig), norm([3*across 2*up]-orig), norm([2*across 2*up]-orig), norm([1*across 2*up]-orig), norm([0 2*up]-orig);...
           norm([5*across 3*up]-orig), norm([4*across 3*up]-orig), norm([3*across 3*up]-orig), norm([2*across 3*up]-orig), norm([1*across 3*up]-orig), norm([0 3*up]-orig);...
           norm([5*across 4*up]-orig), norm([4*across 4*up]-orig), norm([3*across 4*up]-orig), norm([2*across 4*up]-orig), norm([1*across 4*up]-orig), norm([0 4*up]-orig)];
normdistance1 = distance1 - min(distance1(:))
normdistance1 = normdistance1 ./ max(normdistance1(:)) %normalized distance from 1rd sensor

closeness3=distance1 %opposite of distance measuring 'closeness' closer it is to the sensor gets a larger number
normclose3=normdistance1 %normalized closeness to first sensor
closeness1=distance3 %closeness to 1st sensor
normclose1=normdistance3 %normalized closeness to first sensor

closeness2= [norm([5*across 4*up]-orig), norm([4*across 4*up]-orig), norm([3*across 4*up]-orig), norm([2*across 4*up]-orig), norm([1*across 4*up]-orig), norm([0 4*up]-orig);...
           norm([5*across 3*up]-orig), norm([4*across 3*up]-orig), norm([3*across 3*up]-orig), norm([2*across 3*up]-orig), norm([1*across 3*up]-orig), norm([0 3*up]-orig);...
           norm([5*across 2*up]-orig), norm([4*across 2*up]-orig), norm([3*across 2*up]-orig), norm([2*across 2*up]-orig), norm([1*across 2*up]-orig), norm([0 2*up]-orig);...
           norm([5*across 1*up]-orig), norm([4*across 1*up]-orig), norm([3*across 1*up]-orig), norm([2*across 1*up]-orig), norm([1*across 1*up]-orig), norm([0 1*up]-orig);...
           norm([5*across 0*up]-orig), norm([4*across 0*up]-orig), norm([3*across 0*up]-orig), norm([2*across 0*up]-orig), norm([1*across 0*up]-orig), norm([0 0*up]-orig)];
normclose2 = closeness2 - min(closeness2(:))
normclose2 = normclose2 ./ max(normclose2(:))%normalized closeness to second sensor
%% normal peak values
nolimp = 'data/data_10-19-2020_15-49.mat'; %no limp BQ,CQ,DQ,EQ,ER,ES
slightlimp='data/data_10-31-2020_17-08-30'; %slight limp
severelimp='data/data_10-31-2020_17-22-35'; %severe limp

load(severelimp)
filt_datas = lpf_data(datas);
loc_names = {'A1', 'A5', 'E1'};
Fs = 12800;
filt_datas = lpf_data(datas);
clean_data = clean_envelope(filt_datas,Fs);

impactN = 6; % variable input value

[onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impactN,Fs,loc_names);
% 
% gridclose1 = [normclose1(2,2) normclose1(2,3) normclose1(2,4) normclose1(2,5) normclose1(3,5),normclose1(4,5)]
% gridclose2 = [normclose2(2,2) normclose2(2,3) normclose2(2,4) normclose2(2,5) normclose2(5,3),normclose2(4,5)]
% gridclose3 = [normclose3(2,2) normclose3(2,3) normclose3(2,4) normclose3(2,5) normclose3(5,3),normclose3(4,5)]

gridclose1 = [normdistance1(2,2) normdistance1(2,3) normdistance1(2,4) normdistance1(2,5) normdistance1(3,5),normdistance1(4,5)]
gridclose2 = [normdistance2(2,2) normdistance2(2,3) normdistance2(2,4) normdistance2(2,5) normdistance2(5,3),normdistance2(4,5)]
gridclose3 = [normdistance3(2,2) normdistance3(2,3) normdistance3(2,4) normdistance3(2,5) normdistance3(5,3),normdistance3(4,5)]

normpeak1 = peak_val(1,:).*gridclose1
normpeak2 = peak_val(2,:).*gridclose2
normpeak3 = peak_val(3,:).*gridclose3
normpeak = [normpeak1; normpeak2; normpeak3]
figure
bar(normpeak')
ylabel('normdistance times peak value')
xlabel('step #')
legend('sensor1', 'sensor2', 'sensor3')