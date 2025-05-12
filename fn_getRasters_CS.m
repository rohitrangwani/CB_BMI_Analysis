function fn_getRasters(TimeStamps1,events,revents,Fs,savepath,before_zero,after_zero,tStart,tStop,TimeStamps2,meaOnlyFlag,tetrodeSite)
%----------------------------------------------------------------------------------------------------------------
% Function to get rasters of spiking activity around an event and save it in given directory.
% Author@ Aamir Abbasi
% INPUTS:-
%   -- TimeStamps1: 32x6 cell array containing spike timestamps. If split MEA is implanted, first 16 channels
%                   contains units from Cerebellum and the last 16 channels contains units from M1. Otherwise,
%                   all units are from M1. Spike timestamps are in seconds.
%   -- events: Vector containing timestamps of events like trial start or reach onset in samples.
%   -- revents: Vector containing timestamps of reward onset. Non-rewarded trials are marked with zero. 
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
%----------------------------------------------------------------------------------------------------------------

% Define defaults
if nargin < 6
  before_zero = 4;
  after_zero = 4;
  TimeStamps2 = [];
  meaOnlyFlag = 0;
  tetrodeSite = 0;
  tStart = 1;
  tStop = length(events);
end
if nargin < 8
  TimeStamps2 = [];
  meaOnlyFlag = 0;
  tetrodeSite = 0;
  tStart = 1;
  tStop = length(events);
end
if nargin < 10
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
events = events(events>0); %Remove all trials with no reward
if ~isempty(revents)
    revents = revents(tStart:tStop);
else
    revents  = zeros(1,length(events));
end

% Get rasters for MEA recording
for i = 1:size(TimeStamps1,1)
  for j = 2:size(TimeStamps1,2)
    ts = TimeStamps1{i,j};
    if ~isempty(ts)
      for n = 1:length(events)
        time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
        valid_inds = logical((ts>time(1)).*(ts<time(end)));
        raster = ts(valid_inds)-time(1);
        if revents(n)~=0 
          scatter(gca(h),raster,n*ones(1,length(raster)),30,'g','filled');hold on;
%           scatter(gca(h),revents(n)/Fs-time(1),n*ones(1,length(revents(n))),30,'r','filled','d');
        else
%           scatter(gca(h),raster,n*ones(1,length(raster)),30,[140 140 140]/255,'filled');hold on;
          scatter(gca(h),raster,n*ones(1,length(raster)),40,'|','k','linewidth',1.4);hold on;
        end
      end
      vline(before_zero,'k'); title(['CH-',num2str(i),'-SC-',num2str(j)]);
      if meaOnlyFlag == 1 % 1-16 Ch MEA in Cb 17-32 Ch MEA in M1
        if i<=16
          saveas(gcf,[savepath,'RASTERS_Cb-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        else
          saveas(gcf,[savepath,'RASTERS_M1-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
        end
      else % 32 Ch MEA only in M1
        saveas(gcf,[savepath,'RASTERS_M1-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
      end
    end
  end
end
% 
% % Get rasters for tetrode recording
% if ~isempty(TimeStamps2)
%   %TimeStamps2 = TimeStamps2(1:4:size(TimeStamps2,1),:);
%   for i = 1:size(TimeStamps2,1)
%     for j = 2:size(TimeStamps2,2)
%       ts = TimeStamps2{i,j};
%       if ~isempty(ts)
%         for n = 1:length(events)
%           time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
%           valid_inds = logical((ts>time(1)).*(ts<time(end)));
%           raster = ts(valid_inds)-time(1);
%           if revents(n)~=0
%             scatter(gca(h),raster,n*ones(1,length(raster)),30,'g','filled');hold on;
%             scatter(gca(h),revents(n)/Fs-time(1),n*ones(1,length(revents(n))),30,'r','filled','d');
%           else
%             scatter(gca(h),raster,n*ones(1,length(raster)),30,[140 140 140]/255,'filled');hold on;
%           end
%         end
%         vline(before_zero,'k'); title(['CH-',num2str(i),'-SC-',num2str(j)]);
%         if tetrodeSite == 0
%           saveas(gcf,[savepath,'RASTERS_Cb-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
%         else
%           saveas(gcf,[savepath,'RASTERS_M1-Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
%         end
%       end
%     end
%   end
% end
% close;

