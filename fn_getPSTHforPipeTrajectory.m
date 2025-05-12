function [psth_tp, psth_tn] = fn_getPSTHforPipeTrajectory(TimeStamps,events,Fs,Tp,Tn,tStart,tStop,before_zero,after_zero)
%---------------------------------------------------------------------------------------------------
% Function to get rasters of spiking activity around an event and save it in given directory.
% Author@ Aamir Abbasi
% INPUTS:-
%   -- TimeStamps1: 32x6 cell array containing spike timestamps. If split MEA is implanted, first 16 channels
%                   contains units from Cerebellum and the last 16 channels contains units from M1. Otherwise,
%                   all units are from M1. Spike timestamps are in seconds.
%   -- events: Vector containing timestamps of events like trial start or reach onset in samples.
%   -- revents: Vector containing timestamps of reward onset. Non-rewarded trials are marked with zero.
%   -- Fs: Sampling frequency.
%   -- tStart: Starting trial index (default 1).
%   -- tStop: Last trial index (default length(events)).
%   -- before_zero: Time in seconds to look at spike times before an event (default 4 seconds).
%   -- after_zero:  Time in seconds to look at spike times after an event (default 4 seconds).
%----------------------------------------------------------------------------------------------------

% Define defaults
if nargin < 6
  before_zero = 4;
  after_zero = 4;
  tStart = 1;
  tStop = length(events);
end

% Define valid events
events = events(tStart:tStop);

% Define binning parameters
bin_window = 0.001; % in seconds
box_kernel = ones(1,50);
bins = 0:bin_window:before_zero+after_zero;

% Get Tp Tn timestamps for each sortcode
if ~isnan(Tp)
  TimeStamps_tp = TimeStamps(Tp,:);
  TimeStamps_tp = TimeStamps_tp(~cellfun('isempty',TimeStamps_tp));
end

if ~isnan(Tn)
  TimeStamps_tn = TimeStamps(Tn,:);
  TimeStamps_tn = TimeStamps_tn(~cellfun('isempty',TimeStamps_tn));
end

% Collapse Tp and Tn sortcodes in ascending order in to 1D array
if ~isnan(Tp)
  merged_TimeStamps_tp = [];
  for sc = 1:size(TimeStamps_tp,2)
    merged_TimeStamps_tp = sort([merged_TimeStamps_tp TimeStamps_tp{1,sc}]);
  end
end

if ~isnan(Tn)
  merged_TimeStamps_tn = [];
  for sc = 1:size(TimeStamps_tn,2)
    merged_TimeStamps_tn = sort([merged_TimeStamps_tn TimeStamps_tn{1,sc}]);
  end
end

% Get PSTH on multiunit timestamps for Tp and Tn channels
for n = 1:length(events)
  
  time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
  
  if ~isnan(Tp)
    valid_inds = logical((merged_TimeStamps_tp>time(1)).*(merged_TimeStamps_tp<time(end)));
    raster = merged_TimeStamps_tp(valid_inds)-time(1);
    psth_tp(n,:) = conv(histcounts(raster,bins),box_kernel,'same')*(1000/length(box_kernel));
  else
    psth_tp = [];
  end
  
  if ~isnan(Tn)
    valid_inds = logical((merged_TimeStamps_tn>time(1)).*(merged_TimeStamps_tn<time(end)));
    raster = merged_TimeStamps_tn(valid_inds)-time(1);
    psth_tn(n,:) = conv(histcounts(raster,bins),box_kernel,'same')*(1000/length(box_kernel));
  else
    psth_tn = [];
  end
  
end
