%% THIS SCRIPT GENERATES PSTH AND RASTERS FOR INDIRECT UNITS 
%  Author - Aamir Abbasi
% Modified - Rohit
%  Generates .mat files and .tiff plots
%  -----------------------------------------------------------------------
%% PSTHs of indirect units in M1 (Around task-start)
clear;clc;close all;
disp('running');
rootpath = 'Z:\Rohit\BMI_Data\I107\Data\';

bmiBlocks = dir(rootpath);

% tstart = [ 0, 0, 1, 5,  72, 1];
% tstop  = [ 0, 0, 54, 94, 146, 121];
tstart = [0,0, 1, 6, 14, 5];
tstop  = [0,0, 94, 106, 114, 106];

Fs = 1.017252624511719e+03; 
before = 3; %sec
after  = 3;   %sec
for i=3:length(bmiBlocks)
  
  disp(bmiBlocks(i).name);
  
  % Read data
  load([rootpath, bmiBlocks(i).name,'\Timestamps_M1.mat']);
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
  savepath = [rootpath,bmiBlocks(i).name,'\Figs\',bmiBlocks(i).name,'\Indirect_Units_M1\'];

  % Get only good units
  TimeStamps = cell(size(Labels1,1),size(Labels1,2));
  for j=1:size(Labels1,1)
      for k=2:size(Labels1,2)
         if strcmp(Labels1{j,k},'good') == 1 || strcmp(Labels1{j,k},'mua') == 1
             TimeStamps(j,k) = TimeStamps1(j,k);
         end        
      end
  end

  % PSTH 
  [PSTHindirect,~,BINSindirect] = fn_getPSTH_smooth(TimeStamps,all_trials,Fs,savepath,before,after,tstart(i),tstop(i));
  
  % Rasters 
  fn_getRasters(TimeStamps,all_trials,rewards_onset,Fs,savepath,before,after,tstart(i),tstop(i));
     
  % Save
  save([rootpath,bmiBlocks(i).name,'\Indirect_Units_M1.mat'],'PSTHindirect','BINSindirect');
end
disp('Done!!!');

%% PSTHs of indirect units in M1 (Around reward)
clear;clc;close all;
disp('running');

%load BMI sessions info
bmi_session_info;

Fs = 1.017252624511719e+03;
before = 2; %sec
after = 1;   %sec

for s=[12]%M1%15%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-
  
  disp(bmiBlocks(n).name);
  
  % Read data
  load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
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
  savepath = [rootpath,sub{s},'\Data\', bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Indirect_Units_Reward_M1\'];

  % Get only good units
  TimeStamps = cell(size(Labels1,1),size(Labels1,2));
  for j=1:size(Labels1,1)
      for k=2:size(Labels1,2)
         if strcmp(Labels1{j,k},'good') == 1 || strcmp(Labels1{j,k},'mua') == 1
             TimeStamps(j,k) = TimeStamps1(j,k);
         end
      end
  end

  % PSTH 
  [PSTHindirect,~,BINSindirect] = fn_getPSTH_bar(TimeStamps,rewards_onset,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  % Rasters 
  fn_getRasters(TimeStamps,rewards_onset,[],Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
    
  % Save
%   save([rootpath,bmiBlocks(i).name,'\Indirect_Units_Reward_M1.mat'],'PSTHindirect','BINSindirect');
    end
end

disp('Done!!!');

%%