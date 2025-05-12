function [lenLFP_M1, lenLFP_Cb] = fn_getLFPLength(filepath)
%----------------------------------------------------------------------------
% Function to get the duration of TDT block recording by 
% reading the length of the LFP signals. 
% Author@ Aamir Abbasi
% INPUT -
%   -- filepath: file path of the TDT block where .mat files are located
% OUTPUT - 
%   -- lenLFP_M1: Duration of M1 LEF signals
%   -- lenLFP_M1: Duration of Cb LEF signals 
%----------------------------------------------------------------------------

fileNames = {'LFP_M1.mat','LFP_Cb.mat'};
for j = 1:length(fileNames)
  file = fileNames{j};
  load([filepath,'\',file],'l*','fs');
  if exist('lfp_Cb','var') == 1
    lenLFP_Cb = length(lfp_Cb)/fs;
  end
  if exist('lfp_M1','var') == 1
    lenLFP_M1 = length(lfp_M1)/fs;
  end
end

