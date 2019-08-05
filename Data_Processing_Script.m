% *************************************************************************
% This script will process the s2p files. The issue of a limited number of 
% letter labels has been fixed by using letters function written in the
% letter.m file. Which MUST be in the same folder as this script. There is
% no need to open the letters.m file.
% Last updated 8/1/19
% By Jesus Perez
%**************************************************************************
clear
clc
%**************************************************************************
% these are some constants for the signal properties and for surface
% roughness correction.
length=.036; % length of wave guide in millimeters
uo = pi*4e-7; % permeability of free space (H/m)
ur = 1; %relative permeability of copper conductor
s =  5.99e7; % conductivity of copper (S/m)
g = .02; % loss tangent for FR4
E = 4; % dielectic constant for FR4
%**************************************************************************
% create the excel column labels 
Array1 = 1:4000;
newLabels = letters(Array1);
%**************************************************************************
% ask the user for the folder address, a name for the summary folder and
% the height of the average surface roughness in microns
x = inputdlg({'Folder Address','Summary file name','Roughness height (microns)'},'User inputs', [1 100; 1 50; 1 30]);
folder =  x{1};
outputfilename = x{2};
outputfile = [outputfilename,'.xlsx'];
rms = str2double(x{3}); % average copper roughness
rms = rms * 1e-6; % convert micron to meters
%**************************************************************************
% read the files in the folder 
cd(folder)
fnames=dir([folder,'/*.s2p']);
[m,p]=size(fnames);
%trdata = cell(1,m);

K=1;
    fname=fnames(K).name;
    A=importdata(fname,' ',9);
    data= A.data;
    
dataout(:,1)=data(:,1);
dataout(:,2)=data(:,1)*2*pi;
%**************************************************************************
% write the the frequency and omega columns to all the sheets in the
% summary file
warning( 'off', 'MATLAB:xlswrite:AddSheet' ); % remove warning!
disp('Beginning to create the summary file...')
disp(' ')

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

%create RLGC sheets in the summary folder 
xlswrite(outputfile,{'frequency'},'R','a1');
xlswrite(outputfile,dataout(:,1),'R','a2');
xlswrite(outputfile,{'omega'},'R','b1');
xlswrite(outputfile,dataout(:,2),'R','b2'); 

xlswrite(outputfile,{'frequency'},'L','a1');
xlswrite(outputfile,dataout(:,1),'L','a2');
xlswrite(outputfile,{'omega'},'L','b1');
xlswrite(outputfile,dataout(:,2),'L','b2'); 

xlswrite(outputfile,{'frequency'},'G','a1');
xlswrite(outputfile,dataout(:,1),'G','a2');
xlswrite(outputfile,{'omega'},'G','b1');
xlswrite(outputfile,dataout(:,2),'G','b2'); 

xlswrite(outputfile,{'frequency'},'C','a1');
xlswrite(outputfile,dataout(:,1),'C','a2');
xlswrite(outputfile,{'omega'},'C','b1');
xlswrite(outputfile,dataout(:,2),'C','b2'); 

xlswrite(outputfile,{'frequency'},'CORRECTED_ATT_CONST','a1');
xlswrite(outputfile,dataout(:,1),'CORRECTED_ATT_CONST','a2');
xlswrite(outputfile,{'omega'},'CORRECTED_ATT_CONST','b1');
xlswrite(outputfile,dataout(:,2),'CORRECTED_ATT_CONST','b2');

disp('Labels have been created. About to make the individual files and add to the Summary file.')
disp(' ')
disp('Patience,')
disp(' ')
fprintf('There are %d files to process... \n',round(m))
disp(' ')
%**************************************************************************
% Begin looping through all the s2p files and calculate the electrical
% properties and write them into the summary file and create individual
% .xlsx files for them
for k=1:m

    fname=fnames(k).name;
    A=importdata(fname,' ',9);
    data=A.data;

dataout(:,1)=data(:,1);
freq = dataout(:,1);
dataout(:,2)=data(:,1)*2*pi;
w = data(:,1)*2*pi;
dataout(:,3:10)=data(:,2:9);

% Phase unwrap
S11p=dataout(:,4);
S21p=dataout(:,6);
omega=dataout(:,2);

S11up=unwrap(S11p*pi/180);
S21up=unwrap(S21p*pi/180);
% group delay
S21gd=-diff(S21up)./diff(omega).*1000000;
S11gd=-diff(S11up)./diff(omega).*1000000;
%attenuation constants ------>[]
% Combines magnitude and phase values into complex forms
S11=(dataout(:,3).*(exp(1i.*dataout(:,4)/180.*pi)));
S21=(dataout(:,5).*(exp(1i.*dataout(:,6)/180.*pi)));
S12=(dataout(:,7).*(exp(1i.*dataout(:,8)/180.*pi)));
S22=(dataout(:,9).*(exp(1i.*dataout(:,10)/180.*pi)));
% calculate Delta S and define z0
Delta_S= ((S11.*S22)-(S21.*S12));
z0 = 50;
% convert to abcd form
A = (((1+S11).*(1-S22)) + ((S12).*(S21)))./(2.*S21);
B = (z0.*(((1+S11).*(1+S22)) - ((S12).*(S21))))./(1.*S21);
C = (((1-S11).*(1-S22)) - ((S12).*(S21)))./(z0.*2.*S21);
D = (((1-S11).*(1+S22)) + ((S12).*(S21)))./(2.*S21);

Propagation_constant = asinh((B.*C).^(1./2))./length;
att_constant=real(Propagation_constant);
phase_constant=imag(Propagation_constant);

% correct the attentuation constant with the roughness coefficient
D = 1./sqrt(s*pi*uo*ur*freq); % skin depth
x = (rms./D).^2;
Ksr = 1+(2*atan(1.4*x))./pi; % copper surface roughness correction factor

a_diel = 2.3* freq*tan(g)*sqrt(E);% dielectric attentuation
a_con = att_constant-a_diel; % conductor attentuation without surface roughness correction
a_cr = a_con.*Ksr; % conductor attentuation with surface roughness correction

%******************************The RLGC values*****************************
res_top = ((1+S11).*(1+S22)) - ((S12).*(S21));
res_bottom = ((2).*(S21));
res = (res_top)./(res_bottom);
resistance  = (z0).*(res);
R = real(resistance); % resistance per unit length

induct_top = ((1+S11).*(1+S22)) - ((S12).*(S21));
induct_bottom = (2).*(S21);
induct = (induct_top)./(induct_bottom);
inductance = (z0).*(induct);
L = (imag(inductance))./(w); % inductance per unit length

G_top = (1+S11).*(1-S22) +  ( ((S12).*(S21)) - ((2).*(S21)) );
G_bottom  =  (z0).*(((1+S11).*(1+S22)) - ((S12).*(S21)));
Gtivity = (G_top)./(G_bottom);
G = real(Gtivity); % conductivity per unit length

C_top = ((1+S11).*(1-S22) +  ( ((S12).*(S21)) - ((2).*(S21)) ));
C_bottom  = ((z0).*(((1+S11).*(1+S22)) - ((S12).*(S21))));
Capacity  = (C_top)./(C_bottom);
Cap = (imag(Capacity))./(w); % capacitance per unit length
%**************************************************************************
% write to summary folder
% Excel summary file column names
titles={'frequency (Hz)' 'Omega (w)' 'S11'   'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'  ' '  'R'  'L'  'G'  'C' ' ' 'Corrected Attentuation Constant '};
% lim is the limit of where to print values up to
lim =  numel(S12)+1;
lim = num2str(lim);

xlswrite([fname,'.xlsx'],titles,'Sheet1',['a1:w' lim]);
xlswrite([fname,'.xlsx'],dataout,'Sheet1',['a2:j' lim]);
xlswrite([fname,'.xlsx'],S11up,'Sheet1',['p2:p' lim]);
xlswrite([fname,'.xlsx'],S21up,'Sheet1',['l2:l' lim]);
xlswrite([fname,'.xlsx'],S11gd,'Sheet1',['r2:r' lim]);
xlswrite([fname,'.xlsx'],S21gd,'Sheet1',['n2:n' lim]);
% add RLGC values to individual files
xlswrite([fname,'.xlsx'],R,'Sheet1',['t2:t' lim]);
xlswrite([fname,'.xlsx'],L,'Sheet1',['u2:u' lim]);
xlswrite([fname,'.xlsx'],G,'Sheet1',['v2:v' lim]);
xlswrite([fname,'.xlsx'],Cap,'Sheet1',['w2:w' lim]);
xlswrite([fname,'.xlsx'],a_cr,'Sheet1',['y2:y' lim]);

%writting S21 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,dataout(:,5),'AMP_S21',[newLabels{k+2},'2']);
%writting S21 group delay to summmary table
xlswrite(outputfile,{fname},'GD_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21gd,'GD_S21',[newLabels{k+2},'2']);
%writting S21 unwraped phase to summmary table
xlswrite(outputfile,{fname},'UP_S21',[newLabels{k+2},'1']);
xlswrite(outputfile,S21up,'UP_S21',[newLabels{k+2},'2']);
%writting S11 amplitude to summmary table
xlswrite(outputfile,{fname},'AMP_S11',[newLabels{k+2},'1'])
xlswrite(outputfile,dataout(:,3),'AMP_S11',[newLabels{k+2},'2']);
%writting S11 group delay to summmary table
xlswrite(outputfile,{fname},'GD_S11',[newLabels{k+2},'1']);
xlswrite(outputfile,S11gd,'GD_S11',[newLabels{k+2},'2']);
%writting S11 unwraped phase to summmary table
xlswrite(outputfile,{fname},'UP_S11',[newLabels{k+2},'1']);
xlswrite(outputfile,S11up,'UP_S11',[newLabels{k+2},'2']);
%writting ATT_CONST unwraped phase to summmary table
xlswrite(outputfile,{fname},'ATT_CONST',[newLabels{k+2},'1']);
xlswrite(outputfile,att_constant,'ATT_CONST',[newLabels{k+2},'2']);
%writting Corrected atten_CONST unwraped phase to summmary table
xlswrite(outputfile,{fname},'CORRECTED_ATT_CONST',[newLabels{k+2},'1']);
xlswrite(outputfile,a_cr,'CORRECTED_ATT_CONST',[newLabels{k+2},'2']);
%writting PHASE_CONST unwraped phase to summmary table
xlswrite(outputfile,{fname},'PHASE_CONST',[newLabels{k+2},'1']);
xlswrite(outputfile,phase_constant,'PHASE_CONST',[newLabels{k+2},'2']);
%write the RLGC values
xlswrite(outputfile,{fname},'R',[newLabels{k+2},'1']);
xlswrite(outputfile,R,'R',[newLabels{k+2},'2']);
xlswrite(outputfile,{fname},'L',[newLabels{k+2},'1']);
xlswrite(outputfile,L,'L',[newLabels{k+2},'2']);
xlswrite(outputfile,{fname},'G',[newLabels{k+2},'1']);
xlswrite(outputfile,G,'G',[newLabels{k+2},'2']);
xlswrite(outputfile,{fname},'C',[newLabels{k+2},'1']);
xlswrite(outputfile,Cap,'C',[newLabels{k+2},'2']);

c = k;

        if (k == 1)
            fprintf(' %d sheet has been added to the Summary Folder \n',round(c))
        else
            fprintf(' %d sheets have been added to the Summary Folder \n',round(c))
        end

end

datetime.setDefaultFormats('default','hh:mm a MM/dd/yyyy ')
timeStamp = datetime;
disp(' ')
fprintf('Finished at: %s \n',timeStamp)
disp(' ')


