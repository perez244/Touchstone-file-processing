%***************************************************************************
%This script will create a plot for each constant over each frequency in
%the analysis file. Each line represents a 10 minute time interval. 
%Last Edited: 7/10/19
%By Alexandra Dewey
%***************************************************************************

%---------------------------------------------------------------------------
%Asks the user for the file location and analysis file. [1] 
prompt = 'Enter folder address: ';
dd = inputdlg(prompt,'Folder Location',[1 70]);
dd = cell2mat(dd(1,1));
cd(dd) 
filename = uigetfile('*.xlsx');
%---------------------------------------------------------------------------
%Counts the number of sheets in the workbook and their names
[~,sheetNames] = xlsfinfo(filename);
sheetNum = length(sheetNames);
[rows,col] = size((readmatrix(filename,'Sheet', 2))); %Finds the size of one sheet to create a multidimensional array

dataArr = zeros(rows,col,sheetNum);
%Loads sheets into array
for sheetIt = 2:sheetNum
    dataArr(:,:,sheetIt-1) = readmatrix(filename,'Sheet', sheetIt);
    
end
%Extracts data and frequencies from array
freq = dataArr(1,2:end,end-1);
freq = freq/1e9;
dataArr = dataArr(2:end, 2:end,:);
material = erase(filename, 'Analysis file for ');
material = erase(material, '.xlsx');
%---------------------------------------------------------------------------
%Plots the individual graphs
plotNum = 1;
rows = rows-1;
for plotIt = 1:sheetNum-1
    for plotIt2electricboogaloo = 1:10:rows
        hold on
        legend
        plot(freq,dataArr(plotIt2electricboogaloo,:,plotIt))
        freqName = num2str(freq(plotIt2electricboogaloo));
        paraName = char(sheetNames(plotIt+1));
        paraName = strrep(paraName,'_', ' ');
        plotTitle = [material ' ' paraName ' at ' freqName 'GHz'];
        title(plotTitle);
        xlabel('Time min');
        ylabel(paraName);
        plotNum = plotNum+1;
        
    end
    
end
%---------------------------------------------------------------------------
%Saves the graphs as pngs and matlab files
saveChoice = menu('Save plots?', 'Save', 'Cancel');
switch saveChoice
    case 1
            saveName = erase(material,' ');
            saveName = strcat(material, 'FreqPlot');
            saveas(figure(saveIt), saveName, 'png')
    case 2
end
%close
%clear

%---------------------------------------------------------------------------
%Citations:
%[1] Title: DataAnalysisScript
%    Author: Jesus Perez
%    Date: 2019
