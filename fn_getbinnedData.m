function [data,bins] = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tStart,tStop,bin_window,direct,idx)
% Function to get binned  spiking activity around an event and save it in given directory.
% Author@ Aamir Abbasi
% INPUTS:-
%   -- TimeStamps: Nx6 cell array containing spike timestamps.
%   -- events: Vector containing timestamps of events like trial start in samples.
%   -- Fs: Sampling frequency.
%   -- savepath: Path where the user wants to save the files.
%   -- tStart: Starting trial index (default 1).
%   -- tStop: Last trial index (default length(events)). 
%   -- TimeStamps2: 32x2 cell array containing spike timestamps recorded from the Cerebellum/M1 using tetrodes.
%                   By default this cell array is empty.
%   -- meaOnlyFlag: 0 for recordings where a split MEA was implanted in M1 and Cb (default).
%                   1 for recordings where a 32 Ch MEA was implanted in M1.
%   -- before_zero: Time in seconds to look at spike times before an event (default 4 seconds).
%   -- after_zero:  Time in seconds to look at spike times after an event (default 4 seconds).
% OUTPUTS:- 
%   -- PSTH1: A 32x6 cell array containing PSTH for each units in the TimeStamps (M1) cell array. 
%   -- PSTH2: A 32x6 cell array containing PSTH for each units in the TimeStamps2 (Cb) cell array. 
%   -- bins:  A vector containing bin edges
    % Define valid events 
    
    
    events = events(tStart:tStop);
%     events = events(valid_perf<15); 
%     events = events(events>0); 
%     events = events(events>0);

    bins = 0:bin_window:before_zero+after_zero;

%     disp(direct)
    if direct
        data = [];
        ts = TimeStamps(~cellfun('isempty',TimeStamps));
     
%         size(ts)
        ts = ts(idx);
        ts = cell2mat(ts);
        ts = sort(ts);
        
        if ~isempty(ts)
          length(events)  
          for n = 1:length(events)
            time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
            valid_inds = logical((ts>time(1)).*(ts<time(end)));
            raster = ts(valid_inds)-time(1);
            h(n,:) = histcounts(raster,bins);
          end
%           length(events)
%           n
%           disp("direct")
          data = h; 
        end
    else
%         disp("indirect")
        data = cell(size(TimeStamps,1),size(TimeStamps,2));
%         disp(size(data));
        for i = 1:size(TimeStamps,1)
          for j = 1:size(TimeStamps,2) 
              
%             if(idx(i,j))  %for M1 BMI 
            ts = TimeStamps{i,j};
%             disp(~isempty(ts));
            if ~isempty(ts)
%               disp("Non empty ts");
              
              for n = 1:length(events)
                time = ((events(n)/Fs)-before_zero):((events(n)/Fs)+after_zero);
                valid_inds = logical((ts>time(1)).*(ts<time(end)));
                raster = ts(valid_inds)-time(1);
                h(n,:) = histcounts(raster,bins);
              end
%               disp(h);
              data{i,j} = h; 

            end
%             end %for M1 BMI 
          end
        end
    end
    
end    