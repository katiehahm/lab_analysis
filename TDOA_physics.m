% 6/18/21 trying physics approach to localization

clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\ProcessedData\';
T = readtable([data_root_katie,'06-04_06-10_edited_finaldata_nonan.csv']);
data = table2array(T);
arrival_idx = data(:,1:4);
coords = data(:,13:14);
coords = coords + [3584, 3394]; % making sensor 1 coordinate 0,0)

nrows = length(arrival_idx);
predict_coords = zeros(nrows,2);
L1 = 3394+196;
L2 = 3584+3607;
conversions = zeros(nrows,1);

for i = 1:nrows
    a1 = arrival_idx(i,1);
    a2 = arrival_idx(i,2);
    a3 = arrival_idx(i,3);
    a4 = arrival_idx(i,4);
    
    % testing formulas
%     a1 = 30;
%     a2 = 10*(sqrt(41)-2);
%     a3 = 10*(sqrt(61)-2);
%     a4 = 10*(sqrt(45)-2);
%     L1 = 8;
%     L2 = 10;
    
    d = (a3^2 - a4^2 - a2^2 + a1^2)/(2*(a4+a2-a1-a3));
    ad1 = a1+d;
    ad2 = a2+d;
    ad3 = a3+d;
    ad4 = a4+d;

    quadA = L1*L1*(ad4^2 - ad1^2)^2 + L2*L2*(ad2^2-ad1^2)^2;
    quadB = -2*L1*L1*L2*L2*(ad4^2 + ad2^2);
    quadC = L1*L1*L2^4 + L2*L2*L1^4;
    
    csquared = ( -quadB + sqrt(quadB^2 - 4*quadA*quadC) )/(2*quadA);
    csquared2 = ( -quadB - sqrt(quadB^2 - 4*quadA*quadC) )/(2*quadA);
    % how to choose which one?
    c = sqrt(csquared);
    
    conversions(i) = c;
    x = ( c*c*ad4*ad4 - c*c*ad1*ad1 - L2*L2 )/(-2*L2);
    y = ( c*c*ad2*ad2 - c*c*ad1*ad1 - L1*L1 )/(-2*L1);
%     x = (a4^2 - L2^2 + 2*a4*d - a1^2 - 2*a1*d)/(-2*L2);
%     y = (a2^2 + 2*a2*d - a1^2 - 2*a1*d - L1^2)/(-2*L1);
    predict_coords(i,:) = [x,y];
end

coords = real(coords);
predict_coords = real(predict_coords);

diff_x = coords(:,1) - predict_coords(:,1);
diff_z = coords(:,2) - predict_coords(:,2);
diff_x = diff_x(~isnan(diff_x));
diff_z = diff_z(~isnan(diff_z));
error_x = sqrt(mean( diff_x.^2 ));
error_z = sqrt(mean( diff_z.^2 ));

rmse = sqrt(error_x^2 + error_z^2)

figure;
plot(coords(:,1), predict_coords(:,1),'b.')
title('X coordinate')
xlabel('Target coordinate')
ylabel('Predicted coordinate')

figure;
plot(coords(:,2), predict_coords(:,2),'b.')
title('Z coordinate')
xlabel('Target coordinate')
ylabel('Predicted coordinate')