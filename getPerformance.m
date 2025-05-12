%% THIS SCRIPT IS USED FOR GETTING THE LEARNING CURVES
% Rohit 2021-24
%--------------------------------------------------------
%% Read wave channel TTL pulses and build a performance curve
clc; clear; close all;
rootpath = 'Z:\Rohit\BMI_Data\I172\Data\';
                
bmiBlocks = dir(rootpath);

for i=6%:length(bmiBlocks) %start from 3 to ignore . and .. dirs
  
  disp(bmiBlocks(i).name);
  load([rootpath, bmiBlocks(i).name,'\WAV.mat']);
  if exist('wav','var') 
    WAVE = wav;
    clear wav
  end
  if exist('fs','var') 
    Fs = fs;
    clear fs
  end  
  [performance,all_trials, rtrials_onset,rewards_onset]...
    = fn_getPerformanceTimestamps(WAVE,Fs,0,0);
%   pause;
  save([rootpath,bmiBlocks(i).name,'\Events_Performance_PSTH.mat'],...
    'performance','all_trials','rtrials_onset','rewards_onset');
  close all;
end

%% Classify sessions as Robust Learning (RL) or Poor Learning (PL)
clc; clear; close all;
rootpath = 'Z:\Rohit\BMI_Data\I172\Data\';
                
bmiBlocks = dir(rootpath);
%           

%I096
% tstart = [ 5,  72, 1];
% tstop  = [ 94, 146, 121];
%I107
% tstart = [ 1, 6, 14, 5];
% tstop  = [ 94, 106, 114, 106];
%I110
%tstart = [ 0, 44, 31, 1, 15,35,8,46,13,11,38,1];
%tstop  = [ 0, 149, 91, 100, 59,138,96,139,112,33,96,93];
%I111
% tstart = [ 34, 8,52, 46,21,103,20,26];
% tstop  = [ 133, 68, 125, 86,103,151,262,109];

%I112
% tstart = [ 57, 7, 8, 50,42,55,9,20];
% tstop  = [ 128,123,72, 154,139,111,60,50];

%I122
% tstart = [ 0,12, 9, 5, 11,14];
% tstop  = [ 0,200,120,105, 134,89];

% %I116
% tstart = [ 12, 60, 35, 44,2,9,19];
% tstop  = [ 74,83,84, 80,81,103,85];

%I127
% tstart = [ 0,9,56,10 86, 20,1];
% tstop  = [ 0,109,88,82,145,101,150];

%I117
% tstart = [0,0,1,23,1,24,1,21];
% tstop = [0,0,116,45,43,70,32,82];

%I128
% tstart = [ 20,58,60,1];
% tstop  = [ 80,91,103,107];

%I154
% tstart = [0, 66, 12,39,2, 1, 29,14, 19, 45,17,8,  2, 61, 4, 50, 38,29, 5, 1, 15];
% tstop  = [ 0,109,68,61,10,50,83,127,88,101,55,80,65,117,47,159,76,90,80,35,57 ];

%I160
% tstart = [0, 9, 6, 11, 4, 45,19, 11, 6, 1, 20,43, 17,23];
% tstop  = [0,85,80,93,100,82,140,55,100,38,60,120,60,82];

%I161
% tstart = [0,44, 27, 57, 20, 4, 1,  33, 1, 9,42 ];
% tstop  = [0,62, 83, 75, 48, 61,141,113,35,69,90];

%I170
% tstart = [0,1, 40, 1, 1, 21, 8, 5,  21,7 ];
% tstop  = [0,30,85, 78,75,69,49,80,80,40];

%I172
tstart = [0,10, 1, 17, 1, 1];
tstop  = [0,85,80, 57,100,46];

early = [];
late = [];
for i=5%4:length(bmiBlocks)%5:9%[5,6,8,9,10]%:length(bmiBlocks)]%[4:5,7:length(bmiBlocks)]
  disp(bmiBlocks(i).name);
  load([rootpath, bmiBlocks(i).name,'\Events_Performance_PSTH.mat']);
%   if exist('performance_B1','var')
%      performance = performance_B1;
%      clear performance_B1;
%   elseif exist('performance_B2','var')
%      performance= performance_B2;
%      clear performance_B2;
%   end
  valid_perf = performance(tstart(i-2):tstop(i-2));
%   valid_perf = performance;
  early = valid_perf(1:round(length(valid_perf)/3));
  disp(mean(early));
  late  = valid_perf(round(length(valid_perf)*2/3)+1:end);
  disp(mean(late));
  [h,p_val] = ttest2(early,late,'Tail','right');
  display(p_val);
  mean_perf = [mean(early) mean(late)];
  save([rootpath,bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_perf','h','p_val','mean_perf');
end

%% Generate performance early late bar plots for robust learning sessions
clc;clear;close all;
rootpath = 'Z:\Rohit\BMI_Data\I160\Data\';
                
bmiBlocks = dir(rootpath);

for i=4:length(bmiBlocks)-2
  load([rootpath,bmiBlocks(i+2).name,'\Performance_stats_early_late.mat'],'mean_perf');
  mean_perf_consolidated(i,:) = mean_perf(1:2);
end

figure('Color','white'); hold all;
bar(mean(mean_perf_consolidated));
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)));
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(mean_perf_consolidated','k');
% [~,pVal] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2));

%%

% figure('Color','white'); hold all;
for i=4:length(bmiBlocks)
load([rootpath,bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_perf','mean_perf');
early = valid_perf(1:round(length(valid_perf)/3));
late  = valid_perf(round(length(valid_perf)*2/3)+1:end);
perf = [transpose(early), transpose(late)];
figure('Color','white'); hold all;
bar(mean(perf)); 
errorbar([1,2],mean(perf),std(perf)/sqrt((size(perf,1)-1)));
plot(mean_perf','k');
[h,pVal] = ttest(early,late);


display(pVal)
end
%% Percentage of successful vs unsucessful for robust learning sessions
% clc;clear;
rootpath = 'Z:\Rohit\BMI_Data\I128\Data\';
                
bmiBlocks = dir(rootpath);
            
for i=3:length(bmiBlocks)
  load([rootpath,bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_p*');
  if exist('valid_pref','var') == 1
    valid_perf = valid_pref;
  end
  
  early = valid_perf(1:floor(length(valid_perf)/3));
  late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 
  
  tOutcome(i,1) = sum(early==15)/sum(length(early))*100; % successful early
  tOutcome(i,2) = sum(late==15)/sum(length(late))*100; % successful late
  clear valid_pref valid_perf
end

figure('Color','white'); hold all;
bar(mean(tOutcome));
errorbar([1,2],mean(tOutcome),std(tOutcome)/sqrt((size(tOutcome,1)-1)));
% scatter(ones(size(tOutcome,1),1),tOutcome(:,1));
% scatter(ones(size(tOutcome,1),1)*2,tOutcome(:,2));
plot(tOutcome','k');
ylim([0 70]);
[~,pV] = ttest(tOutcome(:,1),tOutcome(:,2));
display(pV);

%% Generate .tiff files of learning curves for all sessions with highlighted early-late trials
clc;clear;close all;
rootpath = 'Z:\Rohit\BMI_Data\I172\Data\';
                
bmiBlocks = dir(rootpath);
          
% tstart = [ 6,  61, 1];
% tstop  = [ 86, 145, 100];
% tstart = [ 6,  61, 1];
% tstop  = [ 86, 140, 80];

% tstart = [ 5,  72, 1];
% tstop  = [ 94, 146, 121];
%I107
% tstart = [ 1, 6, 14, 5];
% tstop  = [ 94, 106, 114, 106];
%I110
% tstart = [0, 44, 31, 1, 15,35,8,46,13,11,38,1];
% tstop  = [ 0, 149, 91, 100, 59,138,96,139,113,33,96,93];
%I111
% tstart = [ 34, 8,52, 46,21,103,20,26];
% tstop  = [ 133, 68, 125, 86,103,151,262,109];

%I112
% tstart = [ 57, 7, 8, 50,42,55,9,20];
% tstop  = [ 128,123,72, 154,139,111,60,50];

%I122
% tstart = [ 0,12, 9, 5, 11,14];
% tstop  = [ 0,200,120,105, 134,89];

%I127
% tstart = [ 0,9,56,10 86, 20,1];
% tstop  = [ 0,109,88,82,145,101,150];

%I116
% tstart = [ 12, 60, 35, 44,2,9,19];
% tstop  = [ 74,83,84, 80,81,103,85];

%I117
% tstart = [0,1,23,1,24,1,21];
% tstop = [0,116,45,43,70,32,82];

%I128
% tstart = [ 20,58,60,1];
% tstop  = [ 80,91,103,107];


%I154
% tstart = [0, 66, 12,39,2, 1, 29,14, 19, 45,17,8,  2, 61, 4, 50, 38,29, 5, 1, 15];
% tstop  = [ 0,109,68,61,10,50,83,127,88,101,55,80,65,117,47,159,76,90,80,35,57 ];

%I160
% tstart = [0, 9, 6, 11, 4, 45,19, 11, 6, 1, 20,43, 17,23];
% tstop  = [0,85,80,93,100,82,140,55,100,38,60,120,60,82];

%I161
% tstart = [0,44, 27, 57, 20, 4, 1,  33, 1, 9,42 ];
% tstop  = [0,62, 83, 75, 48, 61,141,113,35,69,90];

%I170
% tstart = [0,1, 40, 1, 1, 21, 8, 5,  21,7 ];
% tstop  = [0,30,85, 78,75,69,49,80,80,40];

%I172
tstart = [0,10, 1, 17, 1, 1];
tstop  = [0,85,80, 57,100,46];

 figure('Color','white','Position',[941 414 1084 440]);  

for i=5%[4,6:8]%4:length(bmiBlocks)%[4:5,7:9]%4:6%[4,6,8]%
  
  disp(bmiBlocks(i).name);

  % Load data
  load([rootpath,bmiBlocks(i).name,'\Events_Performance_PSTH.mat'],'performance*');
%   if exist('performance_B1','var') == 1
%     performance = performance_B1;
%   elseif exist('performance_B2','var') == 1
%     performance = performance_B2;
%   end
  
  % Get valid trials and classify them in to early-late
  valid_perf = performance(tstart(i-2):tstop(i-2));
  early = valid_perf(1:floor(length(valid_perf)/3));
  late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 
  
  % Plot!
  plot(valid_perf,'k-','LineWidth',2); hold on;
  vline(0,'r-'); vline(floor(length(valid_perf)/3),'r-'); 
  vline(length(valid_perf)-floor(length(valid_perf)/3)+1,'b-'); vline(length(valid_perf),'b-'); 
  title(['BLOCK:',bmiBlocks(i).name]);
  
  % Save plot as .tiff file
  savepath = [rootpath,bmiBlocks(i).name,'\Figs\',bmiBlocks(i).name,'\Performance\'];
  if ~exist(savepath, 'dir')
    mkdir(savepath);
  end  
  saveas(gcf,[savepath,'Trace.tiff']);
  
  % Reinitialize!
  pause(1);
  clf;
  clear perf*  
end
close;
disp('done!');

%% Generate performance early late bar plots for poor learning sessions
clc;clear;close all;
rootpath = 'Z:\Rohit\BMI_Data\I128\Data\';
                
bmiBlocks = dir(rootpath);
            
for i=3:length(bmiBlocks)
  load([rootpath,bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'mean_perf');
  mean_perf_consolidated(i,:) = mean_perf(1:2);
end

figure('Color','white'); hold all;
bar(mean(mean_perf_consolidated));
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/(size(mean_perf_consolidated,1)-1),'.');
scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(mean_perf_consolidated','k')
[~,pVal] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2));

display(pVal);

%% Percentage of success for early vs late trials in poor learning sessions
clc;clear;close all;
rootpath = 'Z:\Rohit\BMI_Data\I117\Data\';
                
bmiBlocks = dir(rootpath);
            
% % Early late classification based on visual inspection
% early_trials_end = [21, 30, 21,  70,  40 ];%[39,30,63,18,59,24,30,39,47,21];
% late_trials_end =  [97, 60, 110, 113, 88];%[199,89,116,75,108,53,60,68,95,97];

for i=[4,7,9]%4:6%[4,6,8]%6:7%11:length(bmiBlocks)
  load([rootpath,bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_p*');
  if exist('valid_pref','var') == 1
    valid_perf = valid_pref;
  end
  
  early = valid_perf(1:floor(length(valid_perf)/3));
  late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 
  
  tOutcome(i,1) = sum(early==15)/sum(length(early))*100; % successful early
  tOutcome(i,2) = sum(late==15)/sum(length(late))*100; % successful late
  clear valid_pref valid_perf
end

figure('Color','white'); hold all;
bar(mean(tOutcome));
errorbar([1,2],mean(tOutcome),std(tOutcome)/(size(tOutcome,1)-1),'.');
scatter(ones(size(tOutcome,1),1),tOutcome(:,1));
scatter(ones(size(tOutcome,1),1)*2,tOutcome(:,2));
plot(tOutcome','k')
[~,pV] = ttest(tOutcome(:,1),tOutcome(:,2));

display(pV);

%% Check all performance curves of robust sessions
clc;clear;close all;
rootpath = 'Z:\Rohit\BMI_Data\I170\Data\';
% rootpath = 'Z:\Rohit\BMI_Data\M1_Data\I060\';                
bmiBlocks = dir(rootpath);


% figure('Color','w','Position',[379 79 1584 899]); hold all;          

for i=5%[3:5]%4:length(bmiBlocks)%[4,7,9]%4:6%length(bmiBlocks)-1%[4,6,8]%[3:5,8:length(bmiBlocks)]
  load([rootpath,bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_perf');
%   scatter([1:length(valid_perf)], valid_perf); hold on;
  %   subplot(4,1,i-2);
%   x = 1:1/length(valid_perf)*100:100;
%   if(length(x)~=length(valid_perf))
%       x= [x,100];
%   end
%   plot(x,smooth(valid_perf,25),'Linewidth',2);
  plot(smooth(valid_perf,25),'Linewidth',2); hold on;
  scatter(1:length(valid_perf),valid_perf,5,'filled'); hold on;
  title(bmiBlocks(i).name);
  xlim([0 100]);
  ylim([0 16]);
  ylabel('Time(s)','Fontsize',18);
  xlabel('Trial','Fontsize',18);
  set(gca, 'TickDir', 'out')
end
%% Percentage of successful vs unsucessful for robust learning sessions
% clc;clear;
rootpath = 'Z:\Rohit\BMI_Data\I128\Data\';
                
bmiBlocks = dir(rootpath);
            
for i=3:length(bmiBlocks)
  load([rootpath,bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_p*');
  if exist('valid_pref','var') == 1
    valid_perf = valid_pref;
  end
  
  early = valid_perf(1:floor(length(valid_perf)/3));
  late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 
  
  tOutcome(i,1) = sum(early==15)/sum(length(early))*100; % successful early
  tOutcome(i,2) = sum(late==15)/sum(length(late))*100; % successful late
  clear valid_pref valid_perf
end

figure('Color','white'); hold all;
bar(mean(tOutcome));
errorbar([1,2],mean(tOutcome),std(tOutcome)/sqrt((size(tOutcome,1)-1)));
% scatter(ones(size(tOutcome,1),1),tOutcome(:,1));
% scatter(ones(size(tOutcome,1),1)*2,tOutcome(:,2));
plot(tOutcome','k');
ylim([0 70]);
[~,pV] = ttest(tOutcome(:,1),tOutcome(:,2));
display(pV);

%% Plot all robust sessions
clc;clear;

% Load BMI session info
bmi_session_info;

rootpath = 'Z:\Rohit\BMI_Data\';
% rootpath = 'Z:\Rohit\BMI_Data\M1_Data\I096\';                
bmiBlocks = dir(rootpath);

figure('Color','w','Position',[379 79 1584 899]); hold all;          

count = 0;
trials = [];
avg = {};
for s=1:15
for i=intersection{s}%:length(bmiBlocks)%[4,7,9]%4:6%length(bmiBlocks)-1%[4,6,8]%[3:5,8:length(bmiBlocks)]
    count = count +1 ;
    bmiBlocks = dir([rootpath,sub{s},'\Data\']);
  load([rootpath,sub{s},'\Data\',bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_perf');
%   scatter([1:length(valid_perf)], valid_perf); hold on;
  %   subplot(4,1,i-2);
  trials(count) = length(valid_perf);
  x = 1:1/length(valid_perf)*101:101;
  while(length(x)<length(valid_perf))
      x= [x,101];
  end
  plot(x,smooth(valid_perf,30),'k','Linewidth',1);
  avg{end+1} = [x;valid_perf];
%   scatter(1:length(valid_perf),valid_perf); hold on;
  ylim([0 16]);
  xlim([0 100]);
  ylabel('Time(s)','Fontsize',18);
  xlabel('Trial','Fontsize',18);
end
end
%%
%% Plot all robust sessions : interpolate and expolate
clc;clear;

% Load BMI session info
bmi_session_info;

rootpath = 'Z:\Rohit\BMI_Data\';
% rootpath = 'Z:\Rohit\BMI_Data\M1_Data\I096\';                
bmiBlocks = dir(rootpath);

figure('Color','w','Position',[379 79 1584 899]); hold all;          
count = 0;
avg = [];
xCommon = linspace(0, 100, 100);
ySum = zeros(1, 100);
for s=1:15
for i=intersection{s}%:length(bmiBlocks)%[4,7,9]%4:6%length(bmiBlocks)-1%[4,6,8]%[3:5,8:length(bmiBlocks)]
    count = count +1;
    
    bmiBlocks = dir([rootpath,sub{s},'\Data\']);
  load([rootpath,sub{s},'\Data\',bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_perf');
%   scatter([1:length(valid_perf)], valid_perf); hold on;
  %   subplot(4,1,i-2);
%   x = 1:1/length(valid_perf)*100:100;
%   while(length(x)<length(valid_perf))
%       x= [x,100];
%   end



%   x = normalize(valid_perf',"range"); 
  y = interp1(1:length(valid_perf), valid_perf, xCommon,'linear','extrap');
  plot(y,'Linewidth',2);
  ySum = y + ySum;
%   scatter(1:length(valid_perf),valid_perf); hold on;
  ylim([0 16]);
  xlim([0 100]);
  ylabel('Time(s)','Fontsize',18);
  xlabel('Trial','Fontsize',18);
end
end
yAverage = ySum / count;
plot(yAverage,'k','Linewidth',4)
%% Plot all robust sessions
clc;clear;

% Load BMI session info
bmi_session_info;

rootpath = 'Z:\Rohit\BMI_Data\';
% rootpath = 'Z:\Rohit\BMI_Data\M1_Data\I096\';                
bmiBlocks = dir(rootpath);

figure('Color','w','Position',[379 79 1584 899]); hold all;          

count = 0;
avg = {};
for s=intact
for i=robust_session{s}%:length(bmiBlocks)%[4,7,9]%4:6%length(bmiBlocks)-1%[4,6,8]%[3:5,8:length(bmiBlocks)]
    count = count +1 ;
    bmiBlocks = dir([rootpath,sub{s},'\Data\']);
  load([rootpath,sub{s},'\Data\',bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_perf');
%   scatter([1:length(valid_perf)], valid_perf); hold on;
  %   subplot(4,1,i-2);
  if (length(valid_perf)>100)
      figure;
      plot(smooth(valid_perf,30),'Linewidth',1); hold on;
      scatter(1:length(valid_perf),valid_perf, 5, 'filled'); hold on;
      title([sub{s},'\Data\',bmiBlocks(i).name]);
      ylim([0 16]);
      xlim([0 100]);
      ylabel('Time(s)','Fontsize',18);
      xlabel('Trial','Fontsize',18);
      set(gca, 'TickDir', 'out')
  end
end
end


%% Get total number of trial for all sessions
clc;clear;

% Load BMI session info
bmi_session_info;
count = 1;
for s=stroke
for i=intersection{s}
    total_trial(count) = tstop{s}(i) - tstart{s}(i);
    count =  count +1;
end
end
