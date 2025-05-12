%% THIS SCRIPT ANALYSIS PSTHs FOR INDIRECT UNITS IN M1 
%  Author - Aamir Abbasi
%  -----------------------------------------------------------------------
%% Analyze PSTHs of direct Tp units in M1 (Around task-start)
clear;clc;close all;
disp('running');
rootpath = 'Z:\Aamir\BMI\';

% % For Poor Learning Sessions
% bmiBlocks =  { 'I050\Data\I050-191219-105728'...
%               ,'I050\Data\I050-191223-133408'...
%               ,'I060\Data\I060-200311-114150'...
%               ,'I061\Data\I061-200505-131845'...
%               ,'I061\Data\I061-200507-111109'...
%               ,'I064\Data\I064-200702-110148'...
%               ,'I064\Data\I064-200707-125708'...
%               ,'I064\Data\I064-200708-112756'...
%               ,'I076\Data\I076-201204-121406'};

% For Robust Learning Sessions          
bmiBlocks =  { 'I076\Data\I076-201201-120931'...
              ,'I076\Data\I076-201202-112106'...
              ,'I076\Data\I076-201203-113433'...
              ,'I076\Data\I076-201205-112030'};
            
PSTH_early = [];  
PSTH_late  = [];
marker = [];
for i=1:length(bmiBlocks)
  
  disp(bmiBlocks{i});
  
  % Read data
  load([rootpath, bmiBlocks{i},'\Indirect_Units_Cb.mat']);
  
  % Remove all empty cells
  PSTHindirect =  PSTHindirect(~cellfun('isempty',PSTHindirect));
  PSTHindirect = PSTHindirect(:);    
  
  % Iterate over units and store PSTH data points in a matrix
  for unit = 1:length(PSTHindirect)
      
      bfr = PSTHindirect{unit};
      
      % Get M1 early/late trials ersp
      eT_Cb = squeeze(bfr(1:floor(size(bfr,1)/3),:));
      lT_Cb = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));
          
      PSTH_early = [PSTH_early;mean(eT_Cb)];
      PSTH_late  = [PSTH_late;mean(lT_Cb)];      
  end
  
  % Mark which units belong to which sessions
  marker = [marker;repmat(i,length(PSTHindirect),1)];
end

% Plot average
figure('Color','white','Position',[744 558 1028 420]);
plot(mean(zscore(PSTH_early,[],2)))
hold on; plot(mean(zscore(PSTH_late,[],2)))
xlim([1000,6000]);
vline(2050);
xlabel('Time (ms)');
ylabel('Firing Rate (z-scored)');
legend('Early','Late');

disp('Done!!!');

%% Analyze PSTHs of indirect units in M1 (Around reward)
clear;clc;close all;
disp('running');
rootpath = 'Z:\Aamir\BMI\';

% % For Poor Learning Sessions
% bmiBlocks =  { 'I050\Data\I050-191219-105728'...
%               ,'I050\Data\I050-191223-133408'...
%               ,'I060\Data\I060-200311-114150'...
%               ,'I061\Data\I061-200505-131845'...
%               ,'I061\Data\I061-200507-111109'...
%               ,'I064\Data\I064-200702-110148'...
%               ,'I064\Data\I064-200707-125708'...
%               ,'I064\Data\I064-200708-112756'...
%               ,'I076\Data\I076-201204-121406'};

% For Robust Learning Sessions          
bmiBlocks =  { 'I076\Data\I076-201201-120931'...
              ,'I076\Data\I076-201202-112106'...
              ,'I076\Data\I076-201203-113433'...
              ,'I076\Data\I076-201205-112030'};
            
PSTH_early = [];  
PSTH_late  = [];
marker = [];
for i=1:length(bmiBlocks)
  
  disp(bmiBlocks{i});
  
  % Read data
  load([rootpath, bmiBlocks{i},'\Indirect_Units_Reward_Cb.mat']);
  
  % Remove all empty cells
  PSTHindirect =  PSTHindirect(~cellfun('isempty',PSTHindirect));
  PSTHindirect = PSTHindirect(:);    
  
  % Iterate over units and store PSTH data points in a matrix
  for unit = 1:length(PSTHindirect)
      
      bfr = PSTHindirect{unit};
      
      % Get M1 early/late trials 
      eT_Cb = squeeze(bfr(1:floor(size(bfr,1)/3),:));
      lT_Cb = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));
          
      PSTH_early = [PSTH_early;mean(eT_Cb)];
      PSTH_late  = [PSTH_late;mean(lT_Cb)];      
  end
  
  % Mark which units belong to which sessions
  marker = [marker;repmat(i,length(PSTHindirect),1)];
end

% Plot average
figure('Color','white','Position',[744 558 1028 420]);
plot(mean(zscore(PSTH_early,[],2)))
hold on; plot(mean(zscore(PSTH_late,[],2)))
xlim([6000,11000]);
vline(10000);
xlabel('Time (ms)');
ylabel('Firing Rate (z-scored)');
legend('Early','Late');

disp('Done!!!');

%% EXAMPLE M1 UNIT PSTH AROUND TASK-START
clear; clc; close all;
path = 'Z:\Aamir\BMI\I076\Data\I076-201202-112106\';
load([path,'Indirect_Units_M1.mat']);

% Get the example PSTH
ch = 28;
sc = 4;
bfr = PSTHindirect{ch,sc};

% Get M1 early/late trials
eT_Cb = squeeze(bfr(1:floor(size(bfr,1)/3),:));
lT_Cb = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));

% Plot
figure('Color','white','Position',[744 558 613 352]);

plot(mean(eT_Cb)); hold on;
plot(mean(lT_Cb));
% ylim([0 40]);
xlim([1000,6000]);
vline(2050);
xlabel('Time (ms)');
ylabel('Firing Rate (Hz)');
legend('Early','Late');

disp('Done!!!');

%% EXAMPLE M1 UNIT PSTH AROUND REWARD
clear; clc; close all;
path = 'Z:\Aamir\BMI\I076\Data\I076-201202-112106\';
load([path,'Indirect_Units_Reward_M1.mat']);

% Get the example PSTH
ch = 14;
sc = 4;
bfr = PSTHindirect{ch,sc};

% Get M1 early/late trials
eT_Cb = squeeze(bfr(1:floor(size(bfr,1)/3),:));
lT_Cb = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));

% Plot
figure('Color','white','Position',[744 558 613 352]);

plot(mean(eT_Cb)); hold on;
plot(mean(lT_Cb));
xlim([6000,11000]);
ylim([0 70]);
vline(10000);
xlabel('Time (ms)');
ylabel('Firing Rate (Hz)');
legend('Early','Late');

disp('Done!!!');

%%

