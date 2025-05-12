%% Get CMR rejected and ERP subtracted LFP signals across trials
%  Author - Aamir Abbasi
% Modified - Rohit
%% --------------------------------------------------------------
clear;clc;close all;
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\I128\Data\';

bmiBlocks = dir(rootpath);

for i=3:length(bmiBlocks)
  
  disp(bmiBlocks(i).name);
  
  % Read collated LFP matrix
  load([rootpath, bmiBlocks(i).name,'\','Collated_LFP_new1.mat']);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR M1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Clear variables common with M1
%   clear trial_data_n ersp_data trial_data1_erp 
%   
%   Z-scoring
%   for j=1:size(trial_data1,1)
%     tmp = reshape(squeeze(trial_data1(j,:,:)),1,[]);
%     trial_data_n(j,:,:) = (trial_data1(j,:,:)-mean(tmp))./std(tmp);
%   end
%   
%   Median subtraction
%   med = median(trial_data_n,1);
%   trial_data1_cmr = bsxfun(@minus,trial_data_n,med);
%     
%   ERP subtraction (Channel by Channel)
%   for ch=1:size(trial_data1_cmr,1)
%     single_channel_data = squeeze(trial_data1_cmr(ch,:,:));
%     mean_erp = mean(single_channel_data,2);
%     trial_data1_erp(ch,:,:) = bsxfun(@minus,single_channel_data,mean_erp);
%   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR CEREBELLUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Clear variables common with M1
  clear trial_data_n ersp_data trial_data2_erp
  
  % Z-scoring
  for j=1:size(trial_data2,1)
    tmp = reshape(squeeze(trial_data2(j,:,:)),1,[]);
    trial_data_n(j,:,:) = (trial_data2(j,:,:)-mean(tmp))./std(tmp);
  end
  
  % Median subtraction CMR
  med = median(trial_data_n,1);
  trial_data2_cmr = bsxfun(@minus,trial_data_n,med);
  
  % ERP subtraction (Channel by Channel)
  for ch=1:size(trial_data2_cmr,1)
    single_channel_data = squeeze(trial_data2_cmr(ch,:,:));
    mean_erp = mean(single_channel_data,2);
    trial_data2_erp(ch,:,:) = bsxfun(@minus,single_channel_data,mean_erp);
  end
  
  % Save
  save([rootpath, bmiBlocks(i).name,'\','Collated_LFP_CMR_ERP.mat'],...%'trial_data1_cmr','trial_data1_erp',
    'trial_data2_cmr','trial_data2_erp',...
    'valid_performance','valid_rtrials','Fs');   
end
disp('done');