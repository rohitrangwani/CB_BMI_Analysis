function [PSTH,bins] = fn_getPSTH_new(TimeStamps,events,Fs,savepath,tStart,tStop,before_zero,after_zero)
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
%   -- before_zero: Time in seconds to look at spike times before an event (default 4 seconds).
%   -- after_zero:  Time in seconds to look at spike times after an event (default 4 seconds).
% OUTPUTS:-
%   -- PSTH1: A 32x6 cell array containing PSTH for each units in the TimeStamps1 (M1) cell array.
%   -- PSTH2: A 32x6 cell array containing PSTH for each units in the TimeStamps2 (Cb) cell array.
%   -- bins:  A vector containing bin edges
%-----------------------------------------------------------------------------------------------------------------------------

% Define defaults
if nargin < 5
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

% Get PSTH for MEA recording
bin_window = 0.001; % in seconds
box_kernel = ones(1,100);
bins = 0:bin_window:before_zero+after_zero;
PSTH = cell(size(TimeStamps,1),size(TimeStamps,2));
for i = 1:size(TimeStamps,1)
  for j = 1:size(TimeStamps,2)
    ts = TimeStamps{i,j};
    if ~isempty(ts)
      for n = 1:length(events)
        time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
        valid_inds = logical((ts>time(1)).*(ts<time(end)));
        raster = ts(valid_inds)-time(1);
        psth(n,:) = conv(histcounts(raster,bins),box_kernel,'same')*(1000/length(box_kernel));
      end
      PSTH{i} = psth; % returned by the function
      
      % for plotting and saving only
      psth_early = mean(psth(1:floor(size(psth,1)/3),:));
      psth_mid   = mean(psth(floor(size(psth,1)/3)+1:floor(size(psth,1)/3)*2,:));
      psth_late  = mean(psth(floor(size(psth,1)/3)*2+1:end,:));
      plot(bins(1:length(psth_early)),psth_early,'b','LineWidth',1.5); hold on;
      plot(bins(1:length(psth_mid)),psth_mid,'g','LineWidth',1.5);
      plot(bins(1:length(psth_late)),psth_late,'m','LineWidth',1.5);
      vline(before_zero,'k'); %xlim([2 6]);
      title(['Neuron-',num2str(i)]);
      legend({'Early','Mid','Late'},'Location','northeast');
      saveas(gcf,[savepath,'PSTH_Ch_',num2str(i),'_Sc_',num2str(j),'.tiff']);clf;
    end
  end
end
close;