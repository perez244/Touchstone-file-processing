% trans_Amp_S21 = [];
% trans_Amp_S11 = [];
% trans_Up_S21 = [];
% trans_Up_S11 = [];
% trans_GD_S21 = [];
% trans_GD_S11 = [];
% trans_Att_Const = [];
% trans_Phase_Const = [];
% TransposedMatlabData = struct('trans_Amp_S21',trans_Amp_S21,'trans_Amp_S11',trans_Amp_S11,'trans_Up_S21',trans_Up_S21,'trans_Up_S11',trans_Up_S11,'trans_GD_S21',trans_GD_S21,'trans_GD_S11',trans_GD_S11,'trans_Att_Const',trans_Att_Const,'trans_Phase_Const',trans_Phase_Const);

% TransposedMatlabData.trans_Amp_S21 = num2cell(TransposedMatlabData.trans_Amp_S21);
% TransposedMatlabData.trans_Amp_S11 = num2cell(TransposedMatlabData.trans_Amp_S11);
% TransposedMatlabData.trans_Up_S21 = num2cell(TransposedMatlabData.trans_Up_S21);
% TransposedMatlabData.trans_Up_S11 = num2cell(TransposedMatlabData.trans_Up_S11);
% TransposedMatlabData.trans_GD_S21 = num2cell(TransposedMatlabData.trans_GD_S21);
% TransposedMatlabData.trans_GD_S11 = num2cell(TransposedMatlabData.trans_GD_S11);
% TransposedMatlabData.trans_Att_Const = num2cell(TransposedMatlabData.trans_Att_Const);
% TransposedMatlabData.trans_Phase_Const = num2cell(TransposedMatlabData.trans_Phase_Const);
% %********************************************************************


% selectedFrequencies = [];
% idx = 1;
% 
% for R = rowLocation
%     
% selectedFrequencies(1,idx) = Frequencies(R,1);
% 
% idx = idx +1;
% 
% end
%**************************************************************
% F =  length(chosenFrequencies);
% for i = 1:F
%     
%     MatlabMatrix(1,1+i) = {chosenFrequencies(1,i)};
%     
% end
%*************************************************************
%%loop for indexing rows and columns in transposing the structure
%P = length(chosenFrequencies);
% P = length(ColumnTitles);
% for J = 1:P %loop through column
%    
%     C2 = 1;
%     for J2 = rowLocation %loop through row
%        
%        TransposedMatlabData.trans_Amp_S11(J,C2) = cell2mat({MatlabData.Amp_S11(J2,J)});
%         
%         C2 = C2+1;
%     end
%     
% end

% %**************************************************************
%for loop for inputing the "transposed" file names
% T = length(ColumnTitles);
% for t = 1:T
%     
%     MatlabMatrix(t,1) = {ColumnTitles{1,t}};
% 
% end
% %**************************************************************************************
% select paramters which you'd like to see, then loop through and match
% their index to the selected parameter

% FN = fieldnames(MatlabData);
% fn = fieldnames(TransposedMatlabData);
% 
% [parametersIndex,tf2] = listdlg('PromptString','Select paramters','ListString',FN);
% 
% % for loop for matching parameters indices to parameter values
% chosenParameters = cell(1,length(parametersIndex));
% z = 1;
% 
% for  r = parametersIndex
%     
% chosenParameters{1,z} = FN{r,1};
% z = z +1;
% 
% end
% chosenParameters = chosenParameters';
% %*******************************************
% neededTransParameters = cell(1,length(parametersIndex));
% y = 1;
% 
% for x = parametersIndex
%     
%     neededTransParameters{1,y} = fn{x,1};
%     y = y + 1;
%     
% end
% 
% neededTransParameters = neededTransParameters';
% % 
% % %*****************************************************************************
% %%%%loop through fields in a struct
% % FN = fieldnames(MatlabData);
% % for Q=1:numel(FN)
% %     if( isnumeric(MatlabData.(FN{Q})))
% %         % do stuff
% %     end
% % end
% % 
% % %%%%%%
% % 
% 
% % FN = fieldnames(MatlabData);
% % fn = fieldnames(TransposedMatlabData);
% 
% for U=1:numel(neededTransParameters)
%     if( isnumeric(TransposedMatlabData.(neededTransParameters{U})) )
%         P = length(chosenFrequencies);
%         for j = 2:P %loop through column
%             c2 = 2;
%             for j2 = rowLocation %loop through row
%                 TransposedMatlabData.(neededTransParameters{U})(j,c2) = cell2mat({MatlabData.(chosenParameters{U})(j2,j)});
%                 c2 = c2+1;
%             end
%         end
%     end
% end


% %*********************************************************************
% for loop for inputting "transposed" Frequencies

% for u=1:numel(neededTransParameters)
%         
%     f =  length(chosenFrequencies);
%     for count = 1:f
%         TransposedMatlabData.(neededTransParameters{u})(1,1+count) = {chosenFrequencies(1,count)};
%     end
%  
% end

% %*******************************************************************
% % for loop for inputing the "transposed" file names
% for X=1:numel(neededTransParameters)
%    
%     T = length(ColumnTitles);
%     for V = 1:T
%         TransposedMatlabData.(neededTransParameters{X})(V,1) = {ColumnTitles{1,V}};
%     
%     end
%     
% end


% %*******************************************************************

% TransposedMatlabData.trans_Amp_S21 = num2cell(TransposedMatlabData.trans_Amp_S21);
% TransposedMatlabData.trans_Amp_S11 = num2cell(TransposedMatlabData.trans_Amp_S11);
% TransposedMatlabData.trans_Up_S21 = num2cell(TransposedMatlabData.trans_Up_S21);
% TransposedMatlabData.trans_Up_S11 = num2cell(TransposedMatlabData.trans_Up_S11);
% TransposedMatlabData.trans_GD_S21 = num2cell(TransposedMatlabData.trans_GD_S21);
% TransposedMatlabData.trans_GD_S11 = num2cell(TransposedMatlabData.trans_GD_S11);
% TransposedMatlabData.trans_Att_Const = num2cell(TransposedMatlabData.trans_Att_Const);
% TransposedMatlabData.trans_Phase_Const = num2cell(TransposedMatlabData.trans_Phase_Const);


%%% now a loop to write to an excel sheet...


% Array = 1:4000;
% Letters = letters(Array);
% NumberOfFiles = length(ColumnTitles)-1;
% % 
% %for loop, for inputing file names
% 
% % %  THis worked!!
% % % box = cell2mat(TransposedMatlabData.trans_Amp_S21)
% % % xlswrite(outputfile,box,'test_dummy','B1')


% % dropping data
% for q =1:numel(neededTransParameters)
%     
%     xlswrite(outputfile,TransposedMatlabData.(neededTransParameters{q}),chosenParameters{q, 1},'B2')
% 
%     
% end



% %for loop for inputting "transposed" labels
% for w=1:numel(neededTransParameters)
%         
%     xlswrite(outputfile,ColumnTitles,chosenParameters{q,1},'A1')
%  
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% dropping frequencies
% 
% for W=1:numel(neededTransParameters)
%         
%     xlswrite(outputfile,chosenFrequencies,chosenParameters{W,1},'B1')
%  
% end
% %*******************************************************************

% %code to create multiple summary files if there are too many s2p files 
% 
% a = 1620; % number of .s2p files we will be crunching 
% b = 500; % fixed number files per Summary file
% r = rem(a,b);
% E = a - r;
% N = E/b; % Number of summary files that will be created
% 
% %create Summary file names. File 1 through file N
% FileIndex = 1:N;
% fileNames = cell(1,N);
% for idx = 1:N
%     fileNames(1,idx)={['Week 2 Board 8 Test Summary File ',int2str(FileIndex(idx)),'.xlsx']};    
% end
% 
% % convert fileName cell to char to use in xlswrite()----> char(fileNames(1,index))
% %**************************************************************************
% % % input the frequencies and omega into all the summary files
% %%% (for) loop through number of files 
% 
%     if(rem(k,b) == 0)
%         % onto the next?
%         if(jdx == N)
%             %do nothing
%         else
%             %filenameIndex = filenameindex + 1
%             jdx = jdx +1;
%         end
%     end
% %**************************************************************************
% code for xlswriting into the different summary folders
index = 1;
sFile = 1;
for k=1:m

    if(rem(k,b) == 0)% b = 5
        %check
        if(index == N)% N=3
            %do nothing
        else
            %filenameIndex = filenameindex + 1
            index = index +1;
        end
        
    end
    
OutPutFile = char(fileNames(1,index));
    
fname=fnames(k).name;
A=importdata(fname,' ',9);
data=A.data;
dataout(:,1)=data(:,1);% frequency values
%...
%...
%...
title={'frequency (Hz)' 'Omega (w)' 'S11'   'S11 Phase' 'S21'   'S21 Phase' 'S12'   'S12 Phase' 'S22'   'S22 Phase'  ' '    'S21 Unwrapped Phase'   ' '  'S21 Group Delay'  ' '     'S11 Unwrapped Phase' ' ' 'S11 Group Delay'};

xlswrite([fname,'.xlsx'],title,'Sheet1','a1:r1');
xlswrite([fname,'.xlsx'],dataout,'Sheet1','a2:j802');
%...
%...
%...
%writting S21 amplitude to summmary table
xlswrite(OutPutFile,{fname},'AMP_S21',[newLabels{k+2},'1']);
xlswrite(OutPutFile,dataout(:,5),'AMP_S21',[newLabels{k+2},'2']);



end









