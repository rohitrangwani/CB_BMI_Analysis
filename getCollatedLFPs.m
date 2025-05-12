%% Get collated lfp signals after rejecting bad trials and channels
%  Author - Aamir Abbasi
% Modified - Rohit
%% -----------------------------------------------------------------
clear;clc;close all;
disp('running...');
% rootpath = 'Z:\Rohit\BMI_Data\I096\Data\';
rootpath = 'Z:\Rohit\BMI_Data\I170\Data\';
savepath = 'Z:\Rohit\BMI_Data\I170\Data\';                

bmiBlocks = dir(rootpath);
                   
% tStart = [ 5,  72, 1];
% tStop  = [ 94, 146, 121];
% tStart = [ 0, 0, 5,  72, 1];
% tStop  = [ 0, 0, 94, 146, 121];
% tStart = [ 0, 0, 1,  72, 1];
% tStop  = [ 0, 0, 54, 146, 121];
%I107
% tStart = [ 0,0,1, 6, 14, 5];
% tStop  = [ 0,0,94, 106, 114, 106];

%I110
% tStart = [ 0,0, 44, 31, 1, 15,35,8,46];
% tStop  = [ 0,0, 149, 91, 100, 59,138,96,139];

%I122
% tStart = [ 0,0,0,12, 9, 5, 11,14];
% tStop  = [ 0,0,0,200,120,105, 134,95];

%I116
% tStart = [0,0, 12, 60, 35, 44,2,9,19];
% tStop  = [0,0, 74,83,84, 80,81,103,85];

%I117
% tStart = [0,0,0,1,23,1,24,1,21];
% tStop = [0,0,0,116,45,43,70,32,82];

%I154
% tStart = [0,0,0, 66, 12,39,2, 1, 29,14, 19, 45,17,8,  2, 61, 4, 50, 38,29, 5, 1, 15];
% tStop  = [0,0, 0,109,68,61,10,50,83,127,88,101,55,80,65,117,47,159,76,90,80,35,57 ];

%I160
% tStart = [0,0,0, 9, 6, 11, 4, 45,19, 11, 6, 1, 20,43, 17,23];
% tStop  = [0,0,0,85,80,93,100,82,140,55,100,38,60,120,60,82];

%I161
% tStart = [0,0,0,44, 27, 57, 20, 4, 1,  33, 1, 9,42 ];
% tStop  = [0,0,0,62, 83, 75, 48, 61,141,113,35,69,90];

%I170
tStart = [0,0,0,1, 40, 1, 1, 21, 8, 5,  21,7 ];
tStop  = [0,0,0,30,85, 78,75,69,49,80,80,40];

% for i=3:3
for i=[5,9,10]%4:length(bmiBlocks)    %[8,14,15,16,18,19]%7:length(bmiBlocks)    
  disp(bmiBlocks(i).name);
  
  % Read LFP signals, trial markers and performance
  matFiles = dir([rootpath, bmiBlocks(i).name,'\','LFP*.mat']);
  for k=1:length(matFiles)
    load([rootpath, bmiBlocks(i).name,'\',matFiles(k).name]);
  end
  load([rootpath, bmiBlocks(i).name,'\','Events_Performance_PSTH.mat'],'all_trial*','rtrials_onse*','performanc*');
  
  % Reassign common variable names 
%   if exist('all_trials_B1','var') == 1
%     all_trials = all_trials_B1;
%     clear all_trials_B1;
%   elseif exist('all_trials_B2','var') == 1
%     all_trials = all_trials_B2;
%     clear all_trails_B2
%   end  
%   if exist('rtrials_onset_B1','var') == 1
%     rtrials_onset = rtrials_onset_B1;
%     clear rtrials_onset_B1
%   elseif exist('rtrials_onset_B2','var') == 1
%     rtrials_onset = rtrials_onset_B2;
%     clear rtrials_onset_B2
%   end  
%   if exist('performance_B1','var') == 1
%     performance = performance_B1;
%     clear performance_B1
%   elseif exist('performance_B2','var') == 1
%     performance = performance_B2;
%     clear performance_B2
%   end  
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
  
  badtrials1 = [];
  trial_data1=[];
   
%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR M1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Check for bad LFP channel(s) by visual inspection
%   badchans1 = fn_checkLFPChans(LFPs1,1:size(LFPs1,2));
  figure;
  for x=1:32
      plot(LFPs1(x,:)+(x*0.005)); hold on;
      title("M1 channels (1-32)");
  end
  badchans1 = input('Enter M1 bad channels (as an array): ');
  
  % Collate LFP signals around trial start
  trial_data1 = fn_collateTrialData(LFPs1,valid_trials,round(4*Fs),round(4*Fs));
  
  % Remove bad channels
  trial_data1(badchans1,:,:) = [];
  
  % Check for bad trials after removing bad channel(s)
  badtrials1 = fn_visualizeTrialData(trial_data1,[]);
  
  %%
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR CEREBELLUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Check for bad LFP channel(s) by visual inspection
%   badchans2 = fn_checkLFPChans(LFPs2,1:size(LFPs2,2));
  figure;
  for x=1:32
      plot(LFPs2(x,:)+(x*0.005)); hold on;
      title("Cb channels (1-32)");
  end
  figure;
  for x=33:64
      plot(LFPs2(x,:)+(x*0.005)); hold on;
      title("Cb channels (33-64)");
  end 
  badchans2 = input('Enter Cb bad channels (as an array): ');
  
  % Collate LFP signals around trial start
  trial_data2 = fn_collateTrialData(LFPs2,valid_trials,round(4*Fs),round(4*Fs));
  
  % Remove bad channels
  trial_data2(badchans2,:,:) = [];

  
  % Check for bad trials after removing bad channel(s)
  badtrials2 = fn_visualizeTrialData(trial_data2,[]);
  %%
  % Remove bad trials
  all_badtrials = union(badtrials1,badtrials2);
  trial_data1(:,:,all_badtrials)= [];
  trial_data2(:,:,all_badtrials)= [];
  valid_performance(all_badtrials)=[];
  valid_rtrials(all_badtrials)=[];
    
  % Save collated lfps
  save([savepath, bmiBlocks(i).name,'\Collated_LFP_new1.mat'],...
      'trial_data1','trial_data2','valid_performance',...
      'valid_rtrials','Fs','badtrials1','badtrials2','badchans1','badchans2');  
end
disp('done');