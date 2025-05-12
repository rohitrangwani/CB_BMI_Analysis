%% THIS SCRIPT GENERATES PSTH AND RASTERS FOR DIRECT UNITS 
%  Author - Aamir Abbasi
%  Modified - Rohit
%  Generates .mat files and .tiff plots
%  -----------------------------------------------------------------------
%% PSTHs of direct units (Around task-start)
clear;clc;close all;
disp('running');
bmi_session_info;


Fs = 1.017252624511719e+03;          
before = 2; %sec
after  = 4;   %sec

FR = {};

for s = 11
rootpath = ['Z:\Rohit\BMI_Data\',sub{s},'\Data\';];


bmiBlocks = dir(rootpath);


for i=intersection{s}(end)
  
  disp(bmiBlocks(i).name);
% Define save path
% savepath = ['Z:\Rohit\BMI_Data\Results\Direct_Units_TS\',bmiBlocks(i).name,'\'];  
savepath = ['Z:\Rohit\BMI_Data\Results\Direct_Units_TS\'];  
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
  


  % PSTH Tp
  [PSTHdirect_tp,~,BINSdirect] = fn_getPSTH_smooth(TimeStamps_tp,all_trials,Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));
  
  % Rasters Tp
  fn_getRasters(TimeStamps_tp,all_trials,rewards_onset,Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));

  
  % PSTH Tn
  [PSTHdirect_tn,~,~] = fn_getPSTH_smooth(TimeStamps_tn,all_trials,Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));
  
  % Rasters Tn
  fn_getRasters(TimeStamps_tn,all_trials,rewards_onset,Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));  
%   
  % Save
%   save([rootpath,bmiBlocks(i).name,'\Direct_Units.mat'],'PSTHdirect_tp','PSTHdirect_tn','BINSdirect');
end
end
disp('Done!!!');

%% PSTHs of direct units(Around reward)
clear;clc;close all;
disp('running');
bmi_session_info;
tic
for s = 1:15
rootpath = ['Z:\Rohit\BMI_Data\',sub{s},'\Data\';];

bmiBlocks = dir(rootpath);

% tstart{s} = [0,0, 1, 6, 14, 5];
% tstop{s}  = [0,0, 94, 106, 114, 106];

Fs = 1.017252624511719e+03;
before = 3; %sec
after = 3;   %sec
for i=intersection{s}
  
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
  savepath = [rootpath,bmiBlocks(i).name,'\Figs\',bmiBlocks(i).name,'\Direct_Units_Reward_AR\'];

  % PSTH
  [PSTHdirect_tp,~,BINSdirect] = fn_getPSTH_smooth(TimeStamps_tp,rewards_onset,Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));
  close;
  % Rasters
  fn_getRasters(TimeStamps_tp,rewards_onset,[],Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));
  close;
  % PSTH
  [PSTHdirect_tn,~,~] = fn_getPSTH_smooth(TimeStamps_tn,rewards_onset,Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));
  close;
  % Rasters
  fn_getRasters(TimeStamps_tn,rewards_onset,[],Fs,savepath,before,after,tstart{s}(i),tstop{s}(i));  
  close;
  % Save
  save([rootpath,bmiBlocks(i).name,'\Direct_Units_Reward.mat'],'PSTHdirect_tp','PSTHdirect_tn','BINSdirect');

end
end
toc
disp('Done!!!');

%% Quantize FR from early to late (Task-start) - Tp and Tn 

clear;clc;close all;
disp('running');
bmi_session_info;


Fs = 1.017252624511719e+03;          
before = 1; %sec
after  = 4;   %sec

FR_early = [];
FR_late = [];
for s = 1:15
rootpath = ['Z:\Rohit\BMI_Data\',sub{s},'\Data\';];

count = 1;
bmiBlocks = dir(rootpath);

for n=intersection{s}
  
  disp(bmiBlocks(n).name);
% Define save path
% savepath = ['Z:\Rohit\BMI_Data\Results\Direct_Units_TS\',bmiBlocks(i).name,'\'];  
savepath = ['Z:\Rohit\BMI_Data\Results\Direct_Units_TS\'];  
  % Read data
  load([rootpath, bmiBlocks(n).name,'\Timestamps_Direct.mat']);
  load([rootpath, bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');  
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
  


  % PSTH Tp
  [PSTHdirect_tp,~,BINSdirect] = fn_getPSTH_smooth(TimeStamps_tp,all_trials,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  for i=1:size(PSTHdirect_tp,1)
      for j =1:size(PSTHdirect_tp,2)
          if ~isempty(PSTHdirect_tp{i,j})
              psth =  PSTHdirect_tp{i,j};
              FR_early(s,count) = mean(mean(psth(1:floor(size(psth,1)/3),:)));
              FR_late(s,count) = mean(mean(psth(end-floor(size(psth,1)/3):end,:)));
              count =  count + 1;
          end
      end
  end
  
  
    % PSTH Tn
  [PSTHdirect_tn,~,BINSdirect] = fn_getPSTH_smooth(TimeStamps_tn,all_trials,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  for i=1:size(PSTHdirect_tn,1)
      for j =1:size(PSTHdirect_tn,2)
          if ~isempty(PSTHdirect_tn{i,j})
              psth =  PSTHdirect_tn{i,j};
              FR_early(s,count) = mean(mean(psth(1:floor(size(psth,1)/3),:)));
              FR_late(s,count) = mean(mean(psth(end-floor(size(psth,1)/3):end,:)));
              count =  count + 1;
          end
      end
  end
  
end
end

  
disp('Done!!!');

%% Quantize FR from early to late (Task-start)- Tn only

clear;clc;close all;
disp('running');
bmi_session_info;


Fs = 1.017252624511719e+03;          
before = 1; %sec
after  = 4;   %sec

FR_early = zeros(15,15);
FR_late = zeros(15,15);
for s = 1:15
rootpath = ['Z:\Rohit\BMI_Data\',sub{s},'\Data\';];

count = 1;
bmiBlocks = dir(rootpath);

for n=intersection{s}
  
  disp(bmiBlocks(n).name);
% Define save path
% savepath = ['Z:\Rohit\BMI_Data\Results\Direct_Units_TS\',bmiBlocks(i).name,'\'];  
savepath = ['Z:\Rohit\BMI_Data\Results\Direct_Units_TS\'];  
  % Read data
  load([rootpath, bmiBlocks(n).name,'\Timestamps_Direct.mat']);
  load([rootpath, bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');  
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
  


  % PSTH Tp
  [PSTHdirect_tn,~,BINSdirect] = fn_getPSTH_smooth(TimeStamps_tn,all_trials,Fs,savepath,before,after,tstart{s}(n),tstop{s}(n));
  
  for i=1:size(PSTHdirect_tn,1)
      for j =1:size(PSTHdirect_tn,2)
          if ~isempty(PSTHdirect_tn{i,j})
              psth =  PSTHdirect_tn{i,j};
              FR_early(s,count) = mean(mean(psth(1:floor(size(psth,1)/3),:)));
              FR_late(s,count) = mean(mean(psth(end-floor(size(psth,1)/3):end,:)));
              count =  count + 1;
          end
      end
  end
end
end

  
disp('Done!!!');
%%

F_intact_early =  FR_early(intact,:);
F_intact_late =  FR_late(intact,:);
F_stroke_early =  FR_early(stroke,:);
F_stroke_late =  FR_late(stroke,:);

F_intact_early =  F_intact_early(:);
F_intact_late =  F_intact_late(:);
F_stroke_early =  F_stroke_early(:);
F_stroke_late =  F_stroke_late(:);

F_intact_early = F_intact_early(F_intact_early(:)~=0);
F_intact_late = F_intact_late(F_intact_late(:)~=0);
F_stroke_early = F_stroke_early(F_stroke_late(:)~=0);
F_stroke_late = F_stroke_late(F_stroke_late(:)~=0);

mean([F_intact_early,F_intact_late])
sem([F_intact_early,F_intact_late])
[h,p,~,t]= ttest(F_intact_early,F_intact_late,'tail','left')
mean([F_stroke_early,F_stroke_late])
sem([F_stroke_early,F_stroke_late])
[h,p,~,t]= ttest(F_stroke_early,F_stroke_late,'tail','left')

%% Percentage change from early to late

F_intact_early =  FR_early(intact,:);
F_intact_late =  FR_late(intact,:);
F_stroke_early =  FR_early(stroke,:);
F_stroke_late =  FR_late(stroke,:);

F_intact_early =  F_intact_early(:);
F_intact_late =  F_intact_late(:);
F_stroke_early =  F_stroke_early(:);
F_stroke_late =  F_stroke_late(:);

F_intact_early = F_intact_early(F_intact_early(:)~=0);
F_intact_late = F_intact_late(F_intact_late(:)~=0);
F_stroke_early = F_stroke_early(F_stroke_late(:)~=0);
F_stroke_late = F_stroke_late(F_stroke_late(:)~=0);

FR_intact = [F_intact_early,F_intact_late];
FR_stroke = [F_stroke_early,F_stroke_late];

Intact_norm = bsxfun(@rdivide,FR_intact,FR_intact(:,1));
Stroke_norm = bsxfun(@rdivide,FR_stroke,FR_stroke(:,1));
mean(Intact_norm)
sem(Intact_norm)
[h,p,~,t]= ttest(Intact_norm(:,1),Intact_norm(:,2),'tail','left')
mean(Stroke_norm)
sem(Stroke_norm)
[h,p,~,t]= ttest(Stroke_norm(:,1),Stroke_norm(:,2),'tail','left')
%% linear mixed effect model

%% linear mixed effect model analysis - intact
rat_id = [];
for row=1:(size(FR_early,1))
    if( ~isempty(find(row==intact)) )
        temp= FR_early(row,:);
        temp=temp(temp(:)~=0);
        rat_id = [rat_id, repelem(row,length(temp))];
    end
end

%
tbl = [[F_intact_early;F_intact_late],[zeros(length(F_intact_early),1);...
    ones(length(F_intact_late),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'FR','EL','RatID'});
formula = 'FR ~ EL + (EL | RatID)';
lme_R2 = fitlme(tbl,formula)
%% linear mixed effect model analysis - stroke

rat_id = [];
for row=1:(size(FR_early,1))
    if( ~isempty(find(row==stroke)) )
        temp= FR_early(row,:);
        temp=temp(temp(:)~=0);
        rat_id = [rat_id, repelem(row,length(temp))];
    end
end
%
tbl = [[F_stroke_early;F_stroke_late],[zeros(length(F_stroke_early),1);...
    ones(length(F_stroke_late),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'FR','EL','RatID'});
formula = 'FR ~ EL + (EL | RatID)';
lme_R2 = fitlme(tbl,formula)
