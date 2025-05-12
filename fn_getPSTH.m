function [PSTH1,PSTH2,bins] = fn_getPSTH(TimeStamps1,events,Fs,savepath,before_zero,after_zero,tStart,tStop,TimeStamps2,meaOnlyFlag,tetrodeSite)
%-----------------------------------------------------------------------------------------------------------------------------
% Function to get psth of spiking activity around an event and save it in given directory.
% Author@ Aamir Abbasi
% INPUTS:-
%   -- TimeStamps1: 32x6 cell array containing spike timestamps. If split MEA is implanted, first 16 channels
%                   contains units from Cerebellum and the last 16 channels contains units from M1. Otherwise,
%                   all units are from M1. Spike timestamps are in seconds.
%   -- events: Vector containing timestamps of events like trial start in samples.
%   -- Fs: Sampling frequency.
%   -- savepath: Path where the user wants to save the .tiff files.
%   -- tStart: Starting trial index (default 1).
%   -- tStop: Last trial index (default length(events)). 
%   -- TimeStamps2: 32x2 cell array containing spike timestamps recorded from the Cerebellum/M1 using tetrodes.
%                   By default this cell array is empty.
%   -- meaOnlyFlag: 0 for recordings where a split MEA was implanted in M1 and Cb (default).
%                   1 for recordings where a 32 Ch MEA was implanted in M1.
%   -- before_zero: Time in seconds to look at spike times before an event (default 4 seconds).
%   -- after_zero:  Time in seconds to look at spike times after an event (default 4 seconds).
%   -- tetrodeSite: 0 when tetrode was implated in Cb (default).
%                   1 when tetrode was implated in M1.
% OUTPUTS:- 
%   -- PSTH1: A 32x6 cell array containing PSTH for each units in the TimeStamps1 (M1) cell array. 
%   -- PSTH2: A 32x6 cell array containing PSTH for each units in the TimeStamps2 (Cb) cell array. 
%   -- bins:  A vector containing bin edges
%-----------------------------------------------------------------------------------------------------------------------------

% Define defaults
if nargin < 5
  before_zero = 4;
  after_zero = 4;
  TimeStamps2 = [];
  meaOnlyFlag = 0;
  tetrodeSite = 0;
  tStart = 1;
  tStop = length(events);
end
if nargin < 7
  TimeStamps2 = [];
  meaOnlyFlag = 0;
  tetrodeSite = 0;
  tStart = 1;
  tStop = length(events);
end
if nargin < 9
  TimeStamps2 = [];
  meaOnlyFlag = 0;
  tetrodeSite = 0;
end

% Define figure properties and get a figure handle
h = figure('Color','white','Position',get(0,'Screensize'));

% Check if the savepath exist, else make it
if ~exist(savepath, 'dir')
  mkdir(savepath);
end

% Define valid events 
events = events(tStart:tStop);
events = events(events>0);

% Get PSTH for MEA recording
bin_window = 0.001; % in seconds
box_kernel = ones(1,100);
bins = 0:bin_window:before_zero+after_zero;
PSTH1 = cell(size(TimeStamps1,1),size(TimeStamps1,2));
for i = 1:size(TimeStamps1,1)
  for j = 2:size(TimeStamps1,2)
    ts = TimeStamps1{i,j};
    if ~isempty(ts)
      for n = 1:length(events)
        time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
        valid_inds = logical((ts>time(1)).*(ts<time(end)));
        raster = ts(valid_inds)-time(1);
        psth(n,:) = conv(histcounts(raster,bins),box_kernel,'same')*(1000/length(box_kernel)); 
      end
      PSTH1{i,j} = psth; % returned by the function
      
      % for plotting and saving only
      psth_early = mean(psth(1:floor(size(psth,1)/3),:));
      psth_mid   = mean(psth(floor(size(psth,1)/3)+1:floor(size(psth,1)/3)*2,:));
      psth_late  = mean(psth(floor(size(psth,1)/3)*2+1:end,:));
      plot(gca(h),bins(1:length(psth_early)),psth_early,'b','LineWidth',1.5); hold on;
      plot(gca(h),bins(1:length(psth_mid)),psth_mid,'g','LineWidth',1.5);
      plot(gca(h),bins(1:length(psth_late)),psth_late,'m','LineWidth',1.5);
      vline(before_zero,'k'); %xlim([2 6]); 
      title(['CH-',num2str(i),'-SC-',num2str(j)]);
      legend({'Early','Mid','Late'},'Location','northeast');
      if meaOnlyFlag == 1 % 1-16 Ch MEA in Cb 17-32 Ch MEA in M1
        if i<=16
          saveas(gcf,[savepath,'PSTH_Cb-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        else
          saveas(gcf,[savepath,'PSTH_M1-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        end
      else % 32 Ch MEA only in M1
         if i>32 
          saveas(gcf,[savepath,'PSTH_Cb-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        else
          saveas(gcf,[savepath,'PSTH_M1-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        end
      end
    end
  end
end

% Get psth for tetrode recording
PSTH2 = cell(size(TimeStamps1,1),size(TimeStamps1,2));
if ~isempty(TimeStamps2)
  %TimeStamps2 = TimeStamps2(1:4:size(TimeStamps2,1),:);
  for i = 1:size(TimeStamps2,1)
    for j = 2:size(TimeStamps2,2)
      ts = TimeStamps2{i,j};
      if ~isempty(ts)
        for n = 1:length(events)
          time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
          ts = TimeStamps2{i,j};
          valid_inds = logical((ts>time(1)).*(ts<time(end)));
          raster = ts(valid_inds)-time(1);
          psth(n,:) = conv(histcounts(raster,bins),box_kernel,'same')*(1000/length(box_kernel)); 
        end
        PSTH2{i,j} = psth; % returned by the function
        
        % for plotting and saving only
        psth_early = mean(psth(1:floor(size(psth,1)/3),:));
        psth_mid   = mean(psth(floor(size(psth,1)/3)+1:floor(size(psth,1)/3)*2,:));
        psth_late  = mean(psth(floor(size(psth,1)/3)*2+1:end,:));
        plot(bins(1:length(psth_early)),psth_early,'b','LineWidth',1.5); hold on;
        plot(bins(1:length(psth_mid)),psth_mid,'g','LineWidth',1.5);
        plot(bins(1:length(psth_late)),psth_late,'m','LineWidth',1.5);
        vline(before_zero,'k'); %xlim([2 6]);
        title(['CH-',num2str(i),'-SC-',num2str(j)]);
        legend({'Early','Mid','Late'},'Location','northeast');
        if tetrodeSite == 0
          saveas(gcf,[savepath,'PSTH_Cb-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        else
          saveas(gcf,[savepath,'PSTH_M1-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        end
      end
    end
  end
end
close;

