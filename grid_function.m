function [error1, error2, error3, error4] = grid_function(filenames,impacts, bounceN,mass,height)
avgs1=[];
    avgs2=[];
    avgs3=[];
    avgs4=[];
    data_root_anya = '/Users/anyachase 1/Dropbox (MIT)/Analysis/data/'
for kk = 1:numel(filenames)
    load([data_root_anya,filenames{kk}]);
    filt_datas = lpf_data(datas);
    loc_names = {'A1', 'A5', 'E1', 'E5'};
    Fs = 12800;
    filt_datas = lpf_data(datas);
    clean_data = clean_envelope2(filt_datas,Fs);
    [onset_idx, peak_idx, peak_val] = TDOA2(clean_data,impacts(kk),Fs,loc_names);

    for i=1:bounceN:impacts(kk)
        avg1= mean(peak_val(1,i:i+4));
        avgs1=[avgs1 avg1]
        avg2= mean(peak_val(2,i:i+4));
        avgs2=[avgs2 avg2];
        avg3= mean(peak_val(3,i:i+4));
        avgs3=[avgs3 avg3];
        avg4= mean(peak_val(4,i:i+4));
        avgs4=[avgs4 avg4];
    end
end
PE3=mass*9.81*height;
[X Y]=meshgrid(1:6,1:6);
E3_1 = avgs1.*[120.958440561460,125.262282886458,97.6770362704578,72.4740300630495,43.3640586843428,27.7814365324730,24.8651464950826,29.3096939023352,53.1670001644637,49.9570158727239,57.9998823990827];
E3_2 = avgs2.*[24.8174652237584,30.1754815152511,34.8384761588134,60.4102303937086,75.4506248534812,49.3933138293228,95.3033942305579,69.3709192131610,47.7245433915168,46.4406916819704,57.0026561034786];
E3_3 = avgs3.*[71.3837973253745,156.054159233132,180.263477917902,168.413645081482,169.304721535371,231.391641286342,256.545721875587,141.109585271027,216.780831064657,203.793921406008,131.683997557909];
E3_4 = avgs4.*[255.286795750698,217.853348202881,214.190591660345,197.295306471823,156.928421852659,148.010933570535,112.872766514682,166.884129211486,178.267707734169,118.985831194972,53.3495539024008];

error1_1= abs([PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3]-E3_1);
error1=error1_1./[PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3];
matrix1 =[1 1 1 1 1 1;
    1 1 1 error1(9) error1(8) error1(7);
    error1(1) error1(2) error1(3) error1(4) error1(5) error1(6);
    1 1 1 1 1 error1(10); 
    1 1 1 1 1 error1(11);
    1 1 1 1 1 1];
strng1=[1 1 1 1 1 1 1 1 1 error1(9) error1(8) error1(7) error1(1) error1(2) error1(3) error1(4) error1(5) error1(6) 1 1 1 1 1 error1(10) 1 1 1 1 1 error1(11) 1 1 1 1 1 1];
string1 =mat2cell(num2str(strng1'),ones(6*6,1));
figure
image(matrix1,'CDataMapping','scaled')
colorbar
hold on
textColor = 'white';
text(Y(:)-.5,X(:)+.25, string1,'Color', textColor,'HorizontalAlignment','left')
%calculte the grid lines
grid = .5:1:6.5;
grid1 = [grid;grid];
grid2 = repmat([.5;6.5],1,length(grid))
%plot the grid lines
plot(grid1,grid2,'k')
plot(grid2,grid1,'k')

error2_1= abs([PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3]-E3_2)
error2=error2_1./[PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3]
matrix2 =[1 1 1 1 1 1;
    1 1 1 error2(9) error2(8) error2(7);
    error2(1) error2(2) error2(3) error2(4) error2(5) error2(6);
    1 1 1 1 1 error2(10); 
    1 1 1 1 1 error2(11);
    1 1 1 1 1 1];
strng2 =[1 1 1 1 1 1 1 1 1 error2(9) error2(8) error2(7) error2(1) error2(2) error2(3) error2(4) error2(5) error2(6) 1 1 1 1 1 error2(10) 1 1 1 1 1 error2(11) 1 1 1 1 1 1];
string2 =mat2cell(num2str(strng2'),ones(6*6,1))
figure
image(matrix2,'CDataMapping','scaled')
colorbar
hold on
textColor = 'white';
text(Y(:)-.5,X(:)+.25, string2,'Color', textColor,'HorizontalAlignment','left')
%calculte the grid lines
grid = .5:1:6.5;
grid1 = [grid;grid];
grid2 = repmat([.5;6.5],1,length(grid))
%plot the grid lines
plot(grid1,grid2,'k')
plot(grid2,grid1,'k')

error3_1= abs([PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3]-E3_3)
error3=error3_1./[PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3]
matrix3 =[1 1 1 1 1 1;
    1 1 1 error3(9) error3(8) error3(7);
    error3(1) error3(2) error3(3) error3(4) error3(5) error3(6);
    1 1 1 1 1 error3(10); 
    1 1 1 1 1 error3(11);
    1 1 1 1 1 1];
strng3 =[1 1 1 1 1 1 1 1 1 error3(9) error3(8) error3(7) error3(1) error3(2) error3(3) error3(4) error3(5) error3(6) 1 1 1 1 1 error3(10) 1 1 1 1 1 error3(11) 1 1 1 1 1 1];
string3 =mat2cell(num2str(strng3'),ones(6*6,1))
figure
image(matrix3,'CDataMapping','scaled')
colorbar
textColor = 'white';
text(Y(:)-.5,X(:)+.25, string3,'Color', textColor,'HorizontalAlignment','left')
hold on
%calculte the grid lines
grid = .5:1:6.5;
grid1 = [grid;grid];
grid2 = repmat([.5;6.5],1,length(grid))
%plot the grid lines
plot(grid1,grid2,'k')
plot(grid2,grid1,'k')

error4_1= abs([PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3]-E3_4)
error4=error4_1./[PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3 PE3]
matrix4 =[1 1 1 1 1 1;
    1 1 1 error4(9) error4(8) error4(7);
    error4(1) error4(2) error4(3) error4(4) error4(5) error4(6);
    1 1 1 1 1 error4(10); 
    1 1 1 1 1 error4(11);
    1 1 1 1 1 1];
strng =[1 1 1 1 1 1 1 1 1 error4(9) error4(8) error4(7) error4(1) error4(2) error4(3) error4(4) error4(5) error4(6) 1 1 1 1 1 error4(10) 1 1 1 1 1 error4(11) 1 1 1 1 1 1];
string4 =mat2cell(num2str(strng'),ones(6*6,1))
figure
% string = mat2cell(num2str([1:6*6]'),ones(6*6,1));
image(matrix4,'CDataMapping','scaled')
colorbar
hold on
%insert the labels
textColor = 'white';
text(Y(:)-.5,X(:)+.25, string4,'Color', textColor,'HorizontalAlignment','left')
%text(matrix4(:),,string,'HorizontalAlignment','left')
%calculte the grid lines
grid = .5:1:6.5;
grid1 = [grid;grid];
grid2 = repmat([.5;6.5],1,length(grid))
%plot the grid lines
plot(grid1,grid2,'k')
plot(grid2,grid1,'k')
end