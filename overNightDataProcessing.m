function overNightDataProcessing(Address)
%Allows for queueing SparQDataProcessingCode by entering the address and file
%name as inputs rather than answering dialog options so matlab can continue
%to run after the work day. 
%Last Edited: 7/10/19
%By: Alex Dewey
%--------------------------------------------------------------------------
Address = cellstr(Address);
overNightSparq 
%All files must be added to Matlab's path in order for the function to work
%(From the "Current Folder" tab in Matlab) Right click -> Add to path ->
%Selected folders and subfolders
end