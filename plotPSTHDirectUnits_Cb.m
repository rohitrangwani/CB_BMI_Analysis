%% THIS SCRIPT GENERATES PSTH AND RASTERS FOR DIRECT UNITS 
%  Author - Aamir Abbasi
%  Generates .mat files and .tiff plots
%  Modified: Rohit Rangwani
%  -----------------------------------------------------------------------
%% PSTHs of direct units in CB (Around task-start)
clear;clc;close all;
disp('running');
rootpath = 'Z:\Rohit\BMI_Data\I112\Data\';

bmiBlocks = dir(rootpath);

% tstart = [ 0, 0, 1, 5,  72, 1];
% tstop  = [ 0, 0, 54, 94, 146, 121];
%I107
% tstart = [0,0, 1, 6, 14, 5];
% tstop  = [0,0, 94, 106, 114, 106];

%I110
% tstart = [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1];
% tstop  = [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93];

%I111
% tstart = [ 0,0,34, 8,52, 46,21,103,20,26];
% tstop  = [ 0,0,133, 68, 125, 86,103,151,262,109];

%I112
tstart = [ 0,0,57, 7, 8, 50,42,55,9,20];
tstop  = [ 0,0,128,123,72, 154,139,111,60,50];


Fs = 1.017252624511719e+03;          
before = 3; %sec
after  = 3;   %sec
for i=3:length(bmiBlocks)
  
  disp(bmiBlocks(i).name);
  
  % Read data
  load([rootpath, bmiBlocks(i).name,'\Timestamps_Direct.mat']);
  load([rootpath, bmiBlocks(i).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');  
  if exist('all_trials_B1','var')
     all_trials = all_trials_B1;
     clear all_trials_B1;
  end
  if exist('all_trials_B2','var')
     all_trials = all_trials_B2;
     clear all_trials_B2;     
  end
  if exist('rewards_onset_B1','var')
     rewards_onset = rewards_onset_B1;
     clear rewards_onset_B1;
  end
  if exist('rewards_onset_B2','var')
     rewards_onset = rewards_onset_B2;
     clear rewards_onset_B2;     
  end  
  
  % Define save path
  savepath = [rootpath,bmiBlocks(i).name,'\Figs\',bmiBlocks(i).name,'\Direct_Units_Reward_TS\'];

  % PSTH Tp
  [PSTHdirect_tp,~,BINSdirect] = fn_getPSTH_smooth(TimeStamps_tp,all_trials,Fs,savepath,before,after,tstart(i),tstop(i));
  
  % Rasters Tp
  fn_getRasters(TimeStamps_tp,all_trials,rewards_onset,Fs,savepath,before,after,tstart(i),tstop(i));
  
  % PSTH Tn
  [PSTHdirect_tn,~,~] = fn_getPSTH_smooth(TimeStamps_tn,all_trials,Fs,savepath,before,after,tstart(i),tstop(i));
  
  % Rasters Tn
  fn_getRasters(TimeStamps_tn,all_trials,rewards_onset,Fs,savepath,before,after,tstart(i),tstop(i));  
%   
  % Save
  save([rootpath,bmiBlocks(i).name,'\Direct_Units.mat'],'PSTHdirect_tp','PSTHdirect_tn','BINSdirect');
end
disp('Done!!!');

%% PSTHs of direct units in Cb (Around reward)
clear;clc;close all;
disp('running');
% rootpath = 'Z:\Rohit\BMI_Data\I107\Data\';

% bmiBlocks = dir(rootpath);

%I107
% tstart = [0,0, 1, 6, 14, 5];
% tstop  = [0,0, 94, 106, 114, 106];

% tstart = [ 0, 0, 1, 5,  72, 1];
% tstop  = [ 0, 0, 54, 94, 146, 121];
 
%I110
% tstart = [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1];
% tstop  = [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93];

%I111
% tstart = [ 0,0,34, 8,52, 46,21,103,20,26];
% tstop  = [ 0,0,133, 68, 125, 86,103,151,262,109];

%I112
% tstart = [ 0,0,57, 7, 8, 50,42,55,9,20];
% tstop  = [ 0,0,128,123,72, 154,139,111,60,50];

%get all bmi session info
bmi_session_info;

Fs = 1.017252624511719e+03;
before = 2; %sec
after = 1;   %sec


for s=12%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
  
  
  disp(bmiBlocks(n).name);
  
  % Read data
  load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
  load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');  
  if exist('all_trials_B1','var')
     all_trials = all_trials_B1;
     clear all_trials_B1;
  end
  if exist('all_trials_B2','var')
     all_trials = all_trials_B2;
     clear all_trials_B2;     
  end
  if exist('rewards_onset_B1','var')
     rewards_onset = rewards_onset_B1;
     clear rewards_onset_B1;
  end
  if exist('rewards_onset_B2','var')
     rewards_onset = rewards_onset_B2;
     clear rewards_onset_B2;     
  end  
  
  % Define save path
  savepath = [rootpath,sub{s},'\Data\', bmiBlocks(n).name,'\Figs\Direct_Units_Reward_AR\'];

  % PSTH
  [PSTHdirect_tp,~,BINSdirect] = fn_getPSTH_bar(TimeStamps_tp,rewards_onset,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  % Rasters
  fn_getRasters(TimeStamps_tp,rewards_onset,[],Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  % PSTH
  [PSTHdirect_tn,~,~] = fn_getPSTH_bar(TimeStamps_tn,rewards_onset,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  % Rasters
  fn_getRasters(TimeStamps_tn,rewards_onset,[],Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));  
  
  % Save
  % save([rootpath,sub{s},'\Data\', bmiBlocks(n).name,''\Direct_Units_Reward.mat'],'PSTHdirect_tp','PSTHdirect_tn','BINSdirect');
    end
end
disp('Done!!!');
close all;

%%