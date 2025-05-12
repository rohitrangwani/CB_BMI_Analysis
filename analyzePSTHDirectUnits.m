%% THIS SCRIPT ANALYSIS PSTHs FOR DIRECT UNITS 
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
bmiBlocks =  { 'I050\Data\I050-191218-112504'...
              ,'I050\Data\I050-191220-104050'...
              ,'I050\Data\I050-191221-121617'...
              ,'I060\Data\I060-200310-112339'...
              ,'I060\Data\I060-200312-111249'...
              ,'I060\Data\I060-200313-113905'...
              ,'I060\Data\I060-200314-131648'...
              ,'I061\Data\I061-200506-110632'...
              ,'I061\Data\I061-200508-120338'...
              ,'I061\Data\I061-200509-122650'...
              ,'I064\Data\I064-200701-112912'...
              ,'I064\Data\I064-200706-120643'...
              ,'I076\Data\I076-201201-120931'...
              ,'I076\Data\I076-201202-112106'...
              ,'I076\Data\I076-201203-113433'...
              ,'I076\Data\I076-201205-112030'};
            
PSTH_early = [];  
PSTH_late  = [];
marker = [];
for i=1:length(bmiBlocks)
  
  disp(bmiBlocks{i});
  
  % Read data
  load([rootpath, bmiBlocks{i},'\Direct_Units.mat']);
  
  % Remove all empty cells
  PSTHdirect_tp =  PSTHdirect_tp(~cellfun('isempty',PSTHdirect_tp));
  PSTHdirect_tp = PSTHdirect_tp(:);    
  
  % Iterate over units and store PSTH data points in a matrix
  for unit = 1:length(PSTHdirect_tp)
      
      bfr = PSTHdirect_tp{unit};
      
      % Get M1 early/late trials ersp
      eT_M1 = squeeze(bfr(1:floor(size(bfr,1)/3),:));
      lT_M1 = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));
          
      PSTH_early = [PSTH_early;mean(eT_M1)];
      PSTH_late  = [PSTH_late;mean(lT_M1)];      
  end
  
  % Mark which units belong to which sessions
  marker = [marker;repmat(i,length(PSTHdirect_tp),1)];
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

%% Analyze PSTHs of direct Tp units in M1 (Around reward)
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
bmiBlocks =  { 'I050\Data\I050-191218-112504'...
              ,'I050\Data\I050-191220-104050'...
              ,'I050\Data\I050-191221-121617'...
              ,'I060\Data\I060-200310-112339'...
              ,'I060\Data\I060-200312-111249'...
              ,'I060\Data\I060-200313-113905'...
              ,'I060\Data\I060-200314-131648'...
              ,'I061\Data\I061-200506-110632'...
              ,'I061\Data\I061-200508-120338'...
              ,'I061\Data\I061-200509-122650'...
              ,'I064\Data\I064-200701-112912'...
              ,'I064\Data\I064-200706-120643'...
              ,'I076\Data\I076-201201-120931'...
              ,'I076\Data\I076-201202-112106'...
              ,'I076\Data\I076-201203-113433'...
              ,'I076\Data\I076-201205-112030'};
            
PSTH_early = [];  
PSTH_late  = [];
marker = [];
for i=1:length(bmiBlocks)
  
  disp(bmiBlocks{i});
  
  % Read data
  load([rootpath, bmiBlocks{i},'\Direct_Units_Reward.mat']);
  
  % Remove all empty cells
  PSTHdirect_tp =  PSTHdirect_tp(~cellfun('isempty',PSTHdirect_tp));
  PSTHdirect_tp = PSTHdirect_tp(:);    
  
  % Iterate over units and store PSTH data points in a matrix
  for unit = 1:length(PSTHdirect_tp)
      
      bfr = PSTHdirect_tp{unit};
      
      % Get M1 early/late trials ersp
      eT_M1 = squeeze(bfr(1:floor(size(bfr,1)/3),:));
      lT_M1 = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));
          
      PSTH_early = [PSTH_early;mean(eT_M1)];
      PSTH_late  = [PSTH_late;mean(lT_M1)];      
  end
  
  % Mark which units belong to which sessions
  marker = [marker;repmat(i,length(PSTHdirect_tp),1)];
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

%% Analyze PSTHs of direct Tn units in M1 (Around task-start)
clear;clc;close all;
disp('running');
rootpath = 'Z:\Aamir\BMI\';

% % For Poor Learning Sessions
% bmiBlocks =  { 'I060\Data\I060-200311-114150'...
%               ,'I061\Data\I061-200505-131845'...
%               ,'I061\Data\I061-200507-111109'...
%               ,'I064\Data\I064-200707-125708'};

% For Robust Learning Sessions          
bmiBlocks =  { 'I060\Data\I060-200310-112339'...
              ,'I060\Data\I060-200312-111249'...
              ,'I060\Data\I060-200313-113905'...
              ,'I060\Data\I060-200314-131648'...
              ,'I061\Data\I061-200506-110632'...
              ,'I061\Data\I061-200508-120338'...
              ,'I061\Data\I061-200509-122650'...
              ,'I076\Data\I076-201201-120931'...
              ,'I076\Data\I076-201202-112106'};
         
PSTH_early = [];  
PSTH_late  = [];
marker = [];
for i=1:length(bmiBlocks)
  
  disp(bmiBlocks{i});
  
  % Read data
  load([rootpath, bmiBlocks{i},'\Direct_Units.mat']);
  
  % Remove all empty cells
  PSTHdirect_tn =  PSTHdirect_tn(~cellfun('isempty',PSTHdirect_tn));
  PSTHdirect_tn = PSTHdirect_tn(:);    
  
  % Iterate over units and store PSTH data points in a matrix
  for unit = 1:length(PSTHdirect_tn)
      
      bfr = PSTHdirect_tn{unit};
      
      % Get M1 early/late trials ersp
      eT_M1 = squeeze(bfr(1:floor(size(bfr,1)/3),:));
      lT_M1 = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));
          
      PSTH_early = [PSTH_early;mean(eT_M1)];
      PSTH_late  = [PSTH_late;mean(lT_M1)];      
  end
  
  % Mark which units belong to which sessions
  marker = [marker;repmat(i,length(PSTHdirect_tn),1)];
end

% Plot average
figure('Color','white','Position',[744 558 1028 420]);
plot(mean(zscore(PSTH_early([10:23],:),[],2)));
hold on; plot(mean(zscore(PSTH_late([10:23],:),[],2)));
xlim([1000,6000]);
vline(2050);
xlabel('Time (ms)');
ylabel('Firing Rate (z-scored)');
legend('Early','Late');

disp('Done!!!');

%% Analyze PSTHs of direct Tn units in M1 (Around reward)
clear;clc;close all;
disp('running');
rootpath = 'Z:\Aamir\BMI\';

% % For Poor Learning Sessions
% bmiBlocks =  { 'I060\Data\I060-200311-114150'...
%               ,'I061\Data\I061-200505-131845'...
%               ,'I061\Data\I061-200507-111109'...
%               ,'I064\Data\I064-200707-125708'};
         
% For Robust Learning Sessions          
bmiBlocks =  { 'I060\Data\I060-200310-112339'...
              ,'I060\Data\I060-200312-111249'...
              ,'I060\Data\I060-200313-113905'...
              ,'I060\Data\I060-200314-131648'...
              ,'I061\Data\I061-200506-110632'...
              ,'I061\Data\I061-200508-120338'...
              ,'I061\Data\I061-200509-122650'...
              ,'I076\Data\I076-201201-120931'...
              ,'I076\Data\I076-201202-112106'};          
            
PSTH_early = [];  
PSTH_late  = [];
marker = [];
for i=1:length(bmiBlocks)
  
  disp(bmiBlocks{i});
  
  % Read data
  load([rootpath, bmiBlocks{i},'\Direct_Units_Reward.mat']);
  
  % Remove all empty cells
  PSTHdirect_tn =  PSTHdirect_tn(~cellfun('isempty',PSTHdirect_tn));
  PSTHdirect_tn = PSTHdirect_tn(:);    
  
  % Iterate over units and store PSTH data points in a matrix
  for unit = 1:length(PSTHdirect_tn)
      
      bfr = PSTHdirect_tn{unit};
      
      % Get M1 early/late trials ersp
      eT_M1 = squeeze(bfr(1:floor(size(bfr,1)/3),:));
      lT_M1 = squeeze(bfr(end-floor(size(bfr,1)/3)+1:end,:));
          
      PSTH_early = [PSTH_early;mean(eT_M1)];
      PSTH_late  = [PSTH_late;mean(lT_M1)];      
  end
  
  % Mark which units belong to which sessions
  marker = [marker;repmat(i,length(PSTHdirect_tn),1)];
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