function fn_getRasters_new(TimeStamps,events,revents,Fs,savepath,tStart,tStop,before_zero,after_zero)
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
%   -- before_zero: Time in seconds to look at spike times before an event (default 4 seconds).
%   -- after_zero:  Time in seconds to look at spike times after an event (default 4 seconds).
%----------------------------------------------------------------------------------------------------------------

% Define defaults
if nargin < 6
  before_zero = 4;
  after_zero = 4;
  tStart = 1;
  tStop = length(events);
end

% Define figure properties and get a figure handle
h = figure('Color','white','Position',get(0,'Screensize'));

% Check if the savepath exist, else make it
if ~exist(savepath, 'dir')
  mkdir(savepath);
end

% Define valid events
events = events(tStart:tStop);
revents = revents(tStart:tStop);

% Get rasters for MEA recording
for i = 1:size(TimeStamps,1)
  for j = 1:size(TimeStamps,2)
    ts = TimeStamps{i,j};
    if ~isempty(ts)
      for n = 1:length(events)
        time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
        valid_inds = logical((ts>time(1)).*(ts<time(end)));
        raster = ts(valid_inds)-time(1);
        if revents(n)~=0
          scatter(gca(h),raster(1:2:end),n*ones(1,length(raster(1:2:end))),30,'g','filled');hold on;
          scatter(gca(h),revents(n)/Fs-time(1),n*ones(1,length(revents(n))),30,'r','filled','d');
        else
          scatter(gca(h),raster(1:2:end),n*ones(1,length(raster(1:2:end))),30,[140 140 140]/255,'filled');hold on;
        end
      end
      vline(before_zero,'k'); title(['Neuron-',num2str(i)]);
      saveas(gcf,[savepath,'RASTERS_Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
    end
  end
end
close;

