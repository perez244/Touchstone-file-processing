% *************************************************************************
% This script will process the s2p files. The issue of a limited number of 
% letter labels has been fixed by using letters function written in the
% letter.m file. Which MUST be in the same folder as the script. There is
% no need to open the letters.m file.

%Edited to run as a function instead of a script - Alex Dewey

% Last updated 7/15/19 by Jesus Perez & Yaw Obeng
%**************************************************************************

% define transmission length or length of the sample
length=.031; % length of wave guide in millimeters

% define folder with the raw data
%folder = 'C:\Users\jap11\Desktop\Jesus Perez Summer 2019\Board-8 Week-2 RF Data'
%folder = 'C:\Users\chuy\Documents\MATLAB\TestData for SparQ Data Processing';

folder = cell2mat(Address(1,1));
outputfile='Sparq Summary File.xlsx';

% read the data in the folder
cd(folder)
fnames=dir([folder,'/*.s2p']);
[m,p]=size(fnames);
trdata=cell(1,m);

% write frequency and omega from the first file to the summary file ----->[]

K=1;
    fname=fnames(K).name;
    A=importdata(fname,' ',9);
    data= A.data;
    
dataout(:,1)=data(:,1);

dataout(:,2)=data(:,1)*2*pi;

xlswrite(outputfile,{'frequency'},'AMP_S21','a1');
xlswrite(outputfile,dataout(:,1),'AMP_S21','a2');
xlswrite(outputfile,{'omega'},'AMP_S21','b1');
xlswrite(outputfile,dataout(:,2),'AMP_S21','b2');

xlswrite(outputfile,{'frequency'},'AMP_S11','a1');
xlswrite(outputfile,dataout(:,1),'AMP_S11','a2');
xlswrite(outputfile,{'omega'},'AMP_S11','b1');
xlswrite(outputfile,dataout(:,2),'AMP_S11','b2');

xlswrite(outputfile,{'frequency'},'UP_S21','a1');
xlswrite(outputfile,dataout(:,1),'UP_S21','a2');
xlswrite(outputfile,{'omega'},'UP_S21','b1');
xlswrite(outputfile,dataout(:,2),'UP_S21','b2');

xlswrite(outputfile,{'frequency'},'UP_S11','a1');
xlswrite(outputfile,dataout(:,1),'UP_S11','a2');
xlswrite(outputfile,{'omega'},'UP_S11','b1');
xlswrite(outputfile,dataout(:,2),'UP_S11','b2');

xlswrite(outputfile,{'frequency'},'GD_S21','a1');
xlswrite(outputfile,dataout(:,1),'GD_S21','a2');
xlswrite(outputfile,{'omega'},'GD_S21','b1');
xlswrite(outputfile,dataout(:,2),'GD_S21','b2');

xlswrite(outputfile,{'frequency'},'GD_S11','a1');
xlswrite(outputfile,dataout(:,1),'GD_S11','a2');
xlswrite(outputfile,{'omega'},'GD_S11','b1');
xlswrite(outputfile,dataout(:,2),'GD_S11','b2');
xlswrite(outputfile,dataout(:,2),'GD_S11','b2');

xlswrite(outputfile,{'frequency'},'ATT_CONST','a1');
xlswrite(outputfile,dataout(:,1),'ATT_CONST','a2');
xlswrite(outputfile,{'omega'},'ATT_CONST','b1');
xlswrite(outputfile,dataout(:,2),'ATT_CONST','b2');

xlswrite(outputfile,{'frequency'},'PHASE_CONST','a1');
xlswrite(outputfile,dataout(:,1),'PHASE_CONST','a2');
xlswrite(outputfile,{'omega'},'PHASE_CONST','b1');
xlswrite(outputfile,dataout(:,2),'PHASE_CONST','b2'); 

nums = 1:4000;
newLabels = letters(nums);

disp('Labels have been created. About to make the individual files and add to the Summary file.')
disp(' ')
disp('Patience...')
disp(' ')
fprintf('We have %f files to process... \n',m)
disp(' ')
%write frequency and omega from the first file to the summary file []-----> 

for k=1:m

    fname=fnames(k).name;

    A_struct=importdata(fname,' ',9);

    data=A_struct.data;

dataout(:,1)=data(:,1);

dataout(:,2)=data(:,1)*2*pi;

%dataout(:,3:10)=data(:,2:9);

% Convert real and imaginary to magnnitude and phase angle

S11mag=( (data(:,2)).^2 + (data(:,3).^2) ).^0.5;
S11phase=atan(data(:,2)./data(:,3));
S12mag=( (data(:,4)).^2 + (data(:,5).^2) ).^0.5;
S12phase=atan(data(:,4)./data(:,5));
S21mag=( (data(:,6)).^2 + (data(:,7).^2) ).^0.5;
S21phase=atan(data(:,6)./data(:,7));
S22mag=( (data(:,8)).^2 + (data(:,9).^2) ).^0.5;
S22phase=atan(data(:,8)./data(:,9));

dataout(:,3) = S11mag;
dataout(:,4) = S11phase;
dataout(:,5) = S12mag;
dataout(:,6) = S12phase;
dataout(:,7) = S21mag;
dataout(:,8) = S21phase;
dataout(:,9) = S22mag;
dataout(:,10) = S22phase;

% Phase unwrap
S11p=S11phase;
S21p=S21phase;
omega=dataout(:,2);

S11up=unwrap(S11p*pi/180);
S21up=unwrap(S21p*pi/180);

% group delay
S21gd=-diff(S21up)./diff(omega).*1000000;
S11gd=-diff(S11up)./diff(omega).*1000000;

%attenuation constants ------>[]

% Combines magnitude and phase values into complex forms

S11=(S11mag.*(exp(1i.*S11phase/180.*pi)));
S21=(S21mag.*(exp(1i.*S21phase/180.*pi)));
S12=(S12mag.*(exp(1i.*S12phase/180.*pi)));
S22=(S22mag.*(exp(1i.*S22phase/180.*pi)));

% calculate Delta S and define z0
Delta_S= ((S11.*S22)-(S21.*S12));
z0 = 50;

% convert to abcd form
A= (((1)+ S11-S22-Delta_S)./(2.*S21));
B= (((1+ S11+S22+Delta_S).*z0)./(2.*S21));
C= (1-S11-S22+Delta_S)./(2.*S21.*z0);
D= (1-S11+S22-Delta_S)./(2.*S21);

Propagation_constant = asinh((B.*C).^(1./2))./length;
att_constant=real(Propagation_constant);
phase_constant=imag(Propagation_constant);


% figure, plot(omega./2/pi,att_constant)
% xlim([4e9 1e10])

% look-up table for output

title={'frequency (Hz)' 'Omega (w)' 'S11'   'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'};

xlswrite([fname,'.xlsx'],title,'Sheet1','a1:r1');
xlswrite([fname,'.xlsx'],dataout,'Sheet1','a2:j802');
xlswrite([fname,'.xlsx'],S11up,'Sheet1','p2:p802');
xlswrite([fname,'.xlsx'],S21up,'Sheet1','l2:l802');
xlswrite([fname,'.xlsx'],S11gd,'Sheet1','r2:r802');
xlswrite([fname,'.xlsx'],S21gd,'Sheet1','n2:n802');

%writting S21 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21mag,'AMP_S21',[newLabels{k+2},'2']);

%writting S21 group delay to summmary table
xlswrite(outputfile,{fname},'GD_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21gd,'GD_S21',[newLabels{k+2},'2']);

%writting S21 unwraped phase to summmary table
xlswrite(outputfile,{fname},'UP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21up,'UP_S21',[newLabels{k+2},'2']);
%writting S11 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S11',[newLabels{k+2},'1'])
xlswrite(outputfile,S11mag,'AMP_S11',[newLabels{k+2},'2']);

%writting S11 group delay to summmary table
xlswrite(outputfile,{fname},'GD_S11',[newLabels{k+2},'1']);
xlswrite(outputfile,S11gd,'GD_S11',[newLabels{k+2},'2']);

%writting S11 unwraped phase to summmary table
xlswrite(outputfile,{fname},'UP_S11',[newLabels{k+2},'1']);
xlswrite(outputfile,S11up,'UP_S11',[newLabels{k+2},'2']);

%writting ATT_CONST unwraped phase to summmary table
xlswrite(outputfile,{fname},'ATT_CONST',[newLabels{k+2},'1']);
xlswrite(outputfile,att_constant,'ATT_CONST',[newLabels{k+2},'2']);

%writting PHASE_CONST unwraped phase to summmary table
xlswrite(outputfile,{fname},'PHASE_CONST',[newLabels{k+2},'1']);
xlswrite(outputfile,phase_constant,'PHASE_CONST',[newLabels{k+2},'2']);

c =fix(k);
fprintf(' %f sheets have been added to the Summary Folder \n',c)

end

disp(' ')
disp('Done!')
disp(' ')
