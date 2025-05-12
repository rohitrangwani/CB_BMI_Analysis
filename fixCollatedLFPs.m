
clear;clc;close all;
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\I096\Data\';
                
bmiBlocks = dir(rootpath);
tStart = [ 0, 0, 5,  72, 1];
tStop  = [ 0, 0, 94, 146, 121];

for i=3:length(bmiBlocks)
    
  disp(bmiBlocks(i).name);
  matFiles = dir([rootpath, bmiBlocks(i).name,'\','LFP*.mat']);
  for k=1:length(matFiles)
    load([rootpath, bmiBlocks(i).name,'\',matFiles(k).name]);
  end
  load([rootpath, bmiBlocks(i).name,'\','Events_Performance_PSTH.mat'],'all_trial*','rtrials_onse*','performanc*');
  
  if exist('fs','var') == 1
    Fs = fs;
    clear fs;
  end     
  if exist('lfp_M1','var') == 1
    LFPs1 = lfp_M1;
    clear lfp_M1;
  end      
  if exist('lfp_Cb','var') == 1
    LFPs2 = lfp_Cb;
    clear lfp_Cb;
  end  
     
  % Get valid trials
  valid_trials      = all_trials(tStart(i):tStop(i));
  valid_rtrials     = rtrials_onset(tStart(i):tStop(i));
  valid_performance = performance(tStart(i):tStop(i));
    
  trial_data1 = fn_collateTrialData(LFPs1,valid_trials,round(4*Fs),round(4*Fs));
  trial_data2 = fn_collateTrialData(LFPs2,valid_trials,round(4*Fs),round(4*Fs));
  
  %Load the bad trials and channels
  load([rootpath, bmiBlocks(i).name,'\Collated_LFP_new.mat'],'badtrials*','badchans*');
 
  % Remove bad channels
  trial_data1(badchans1,:,:) = [];
  trial_data2(badchans2,:,:) = [];
  
  % Remove bad trials
  all_badtrials = union(badtrials1,badtrials2);
  trial_data1(:,:,all_badtrials)= [];
  trial_data2(:,:,all_badtrials)= [];
  valid_performance(all_badtrials)=[];
  valid_rtrials(all_badtrials)=[];
    
  % Save collated lfps
  save([rootpath, bmiBlocks(i).name,'\Collated_LFP_new1.mat'],...
      'trial_data1','trial_data2','valid_performance',...
      'valid_rtrials','Fs','badtrials1','badtrials2','badchans1','badchans2');
end

