%% THIS SCRIPT GENERATES PSTH AND RASTERS FOR INDIRECT UNITS 
%  Author - Rohit
%  Generates .mat files and .tiff plots
%  -----------------------------------------------------------------------
%% PSTHs of indirect units in Cb (Around task-start)
clear;clc;close all;
disp('running');
rootpath = 'Z:\Rohit\BMI_Data\I1\Data\';

bmiBlocks = dir(rootpath);

% tstart = [ 0, 0, 1, 5,  72, 1];
% tstop  = [ 0, 0, 54, 94, 146, 121];
%I107
% tstart = [0,0, 1, 6, 14, 5];
% tstop  = [0,0, 94, 106, 114, 106];                             
%I110
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
  load([rootpath, bmiBlocks(i).name,'\Timestamps_Cb.mat']);
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
  savepath = [rootpath,bmiBlocks(i).name,'\Figs\',bmiBlocks(i).name,'\Indirect_Units_Cb_TS\'];

  % Get only good units
  TimeStamps = cell(size(Labels2,1),size(Labels2,2));
  for j=1:size(Labels2,1)
      for k=2:size(Labels2,2)
         if strcmp(Labels2{j,k},'good') == 1
             TimeStamps(j,k) = TimeStamps2(j,k);
         end
      end
  end


  % PSTH 
  [PSTHindirect,~,BINSindirect] = fn_getPSTH_smooth(TimeStamps,rewards_onset,Fs,savepath,before,after,tstart(i),tstop(i));
  [PSTHindirect,~,BINSindirect] = fn_getPSTH(TimeStamps,rewards_onset,Fs,savepath,before,after,tstart(i),tstop(i));
  % Rasters 
  fn_getRasters_sample(TimeStamps,rewards_onset,[],Fs,savepath,before,after,tstart(i),tstop(i));
  fn_getRasters(TimeStamps,rewards_onset,[],Fs,savepath,before,after,tstart(i),tstop(i));
     
  % Save
  save([rootpath,bmiBlocks(i).name,'\Indirect_Units_Cb.mat'],'PSTHindirect','BINSindirect');
end
disp('Done!!!');

%% PSTHs of indirect units in Cb (Around reward)
clear;clc;close all;

Fs = 1.017252624511719e+03;
before = 2.2; %sec
after = 0.2;   %sec
tstart = {[0,0, 1, 6, 14, 5], [ 0, 0, 1, 6,  13, 5],...
    [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1],...
    [ 0,0,34, 8,52, 46,21,103,20,26],...
    [ 0,0,57, 7, 8, 50,42,55,9,20]...
    [ 0,0,0,12, 9, 5, 11,14]...
    [ 0,0,12, 60, 35, 44,2,9,19]...
    [0,0,0,0,1,23,1,24,1,21]...
    [0,0, 0,9,56,10 86, 20,1]...
    [0,0, 20,58,60,1]};

tstop  = {[0,0, 94, 106, 114, 106],[ 0, 0, 104, 105, 114, 105],...
    [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93],...
    [ 0,0,133, 68, 125, 86,103,151,262,109],...
    [ 0,0,128,123,72, 154,139,111,60,50]...
    [0,0, 0,200,120,105, 134,89]...
    [0,0, 74,83,84, 80,81,103,85]...
    [0,0,0,0,116,45,43,70,32,82]...
    [0,0, 0,109,88,82,145,101,150]...
    [0,0, 80,91,103,107]};



sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};


disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6]};

for s=3:10%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1

  
  disp(bmiBlocks(n).name);
  
  % Read data

  load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
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
  savepath = [rootpath, sub{s},'\Data\',bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Indirect_Units_Reward_Cb_AR\'];

  % Get only good units
  TimeStamps = cell(size(Labels2,1),size(Labels2,2));
  for j=1:size(Labels2,1)
      for k=2:size(Labels2,2)
          
         if strcmp(Labels2{j,k},'good') == 1
             TimeStamps(j,k) = TimeStamps2(j,k);
         end
      end
  end

  % PSTH 
%   [PSTHindirect,~,BINSindirect] = fn_getPSTH_bar(TimeStamps,rewards_onset,Fs,savepath,before,after,tstart(i),tstop(i));
  [PSTHindirect,~,BINSindirect] = fn_getPSTH(TimeStamps,rewards_onset,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  % Rasters 
%   fn_getRasters(TimeStamps,rewards_onset,[],Fs,savepath,before,after,tstart(i),tstop(i));
%   fn_getRasters(TimeStamps,rewards_onset,[],Fs,savepath,before,after,tstart(i),tstop(i));
  % Save
  save([rootpath, sub{s},'\Data\',bmiBlocks(n).name,'\Indirect_Units_Reward_Cb_new.mat'],'PSTHindirect','BINSindirect');
    end
end
disp('Done!!!');

%%

%% PSTHs of direct units in Cb (Around reward)
clear;clc;close all;
clear;clc;close all;

Fs = 1.017252624511719e+03;
before = 2.2; %sec
after = 0.2;   %sec
tstart = {[0,0, 1, 6, 14, 5], [ 0, 0, 1, 6,  13, 5],...
    [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1],...
    [ 0,0,34, 8,52, 46,21,103,20,26],...
    [ 0,0,57, 7, 8, 50,42,55,9,20]...
    [ 0,0,0,12, 9, 5, 11,14]...
    [ 0,0,12, 60, 35, 44,2,9,19]...
    [0,0,0,0,1,23,1,24,1,21]...
    [0,0, 0,9,56,10 86, 20,1]...
    [0,0, 20,58,60,1]};

tstop  = {[0,0, 94, 106, 114, 106],[ 0, 0, 104, 105, 114, 105],...
    [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93],...
    [ 0,0,133, 68, 125, 86,103,151,262,109],...
    [ 0,0,128,123,72, 154,139,111,60,50]...
    [0,0, 0,200,120,105, 134,89]...
    [0,0, 74,83,84, 80,81,103,85]...
    [0,0,0,0,116,45,43,70,32,82]...
    [0,0, 0,109,88,82,145,101,150]...
    [0,0, 80,91,103,107]};



sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};


disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6]};

for s=2:10%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1

  
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
  savepath =  [rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];

  % PSTH
  [PSTHdirect_tp,~,BINSdirect] = fn_getPSTH(TimeStamps_tp,rewards_onset,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  % Rasters
%   fn_getRasters(TimeStamps_tp,rewards_onset,[],Fs,savepath,before,after,tstart(i),tstop(i));
  
  % PSTH
  [PSTHdirect_tn,~,~] = fn_getPSTH(TimeStamps_tn,rewards_onset,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  % Rasters
%   fn_getRasters(TimeStamps_tn,rewards_onset,[],Fs,savepath,before,after,tstart(i),tstop(i));  
  
  % Save
  save([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Direct_Units_Reward_new.mat'],'PSTHdirect_tp','PSTHdirect_tn','BINSdirect');
    end
end
disp('Done!!!');