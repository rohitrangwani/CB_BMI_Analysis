function [PSTH1,bins] = fn_getMUPSTH(TimeStamps1,events,Fs,savepath,before_zero,after_zero,tStart,tStop,TpTnflag)
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
%   -- before_zero: Time in seconds to look at spike times before an event (default 4 seconds).
%   -- after_zero:  Time in seconds to look at spike times after an event (default 4 seconds).
%   -- tStart: Starting trial index (default 1).
%   -- tStop: Last trial index (default length(events)).
%   -- TpTnflag: For saving purpose,0:tp 1:tn (default 0).
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
    TpTnflag = 0;
end
if nargin < 7
    tStart = 1;
    tStop = length(events);
    TpTnflag = 0;
end
if nargin < 9
   TpTnflag = 0; 
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
TimeStamps1 = TimeStamps1(~cellfun('isempty',TimeStamps1));
ts = [];
for i = 1:size(TimeStamps1,1)
    ts = [ts;TimeStamps1{i}];
    ts = sort(ts);
end
for n = 1:length(events)
    time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
    valid_inds = logical((ts>time(1)).*(ts<time(end)));
    raster = ts(valid_inds)-time(1);
    psth(n,:) = conv(histcounts(raster,bins),box_kernel,'same')*(1000/length(box_kernel));
end
PSTH1 = psth; % returned by the function

% for plotting and saving only
psth_early = mean(psth(1:floor(size(psth,1)/3),:));
psth_mid   = mean(psth(floor(size(psth,1)/3)+1:floor(size(psth,1)/3)*2,:));
psth_late  = mean(psth(floor(size(psth,1)/3)*2+1:end,:));
plot(gca(h),bins(1:length(psth_early)),psth_early,'b','LineWidth',1.5); hold on;
plot(gca(h),bins(1:length(psth_mid)),psth_mid,'g','LineWidth',1.5);
plot(gca(h),bins(1:length(psth_late)),psth_late,'m','LineWidth',1.5);
vline(before_zero,'k');
% xlim([0 8]);
legend({'Early','Mid','Late'},'Location','northeast');
if TpTnflag == 0
    saveas(gcf,[savepath,'PSTH_M1-Tp-Allunits.tiff']);clf;
end
if TpTnflag == 1
    saveas(gcf,[savepath,'PSTH_M1-Tn-Allunits.tiff']);clf;
end
close;