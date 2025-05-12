%% Data high GPFA analysis - Cb
% Rohit

%%

clc;clear;



sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};

stroke= [4:5,9:10];

intact = [1:10];
intact = setdiff(intact, stroke);

Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';

tstart = {[0,0, 1, 6, 14, 5], [ 0, 0, 1, 6,  13, 5],...
    [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1],...
    [ 0,0,34, 8,52, 46,21,103,20,26],...
    [ 0,0,57, 7, 8, 50,42,55,9,20]...
    [ 0,0,0,12, 9, 5, 11,14]...
    [ 0,0,12, 60, 35, 44,2,9,19]...
    [0,0,0,0,1,23,1,24,1,21]...
    [0,0, 0,9,56,10 86, 20,1]...
    [0,0, 20,58,60,1]};

tstop  = {[0,0, 94, 106, 114, 106],[ 0, 0, 104, 105, 114, 105],...
    [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93],...
    [ 0,0,133, 68, 125, 86,103,151,262,109],...
    [ 0,0,128,123,72, 154,139,111,60,50]...
    [0,0, 0,200,120,105, 134,89]...
    [0,0, 74,83,84, 80,81,103,85]...
    [0,0,0,0,116,45,43,70,32,82]...
    [0,0, 0,109,88,82,145,101,150]...
    [0,0, 80,91,103,107]};

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6]};

% bmiBlocks = dir(rootpath); 
before_zero = 2.2;
after_zero =0.2;

tp = [];
tn = [];

res = [];
indirect = [];

% bin = [5,10,20,25,50,100]; % binnning at different ms
% for b=4%1:length(bin)
for s=1%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=4%robust_session{s}
    
        
    % Load timestamps
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
        bin_window = 25/1000; %in sec
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
        events = rewards_onset;
        valid_perf = performance(tstart{s}(n):tstop{s}(n)); 
                
        TimeStamps = cell(size(Labels2,1),size(Labels2,2));
        
        % get good trials only
        for j=1:size(Labels2,1)
          for k=2:size(Labels2,2)
             if strcmp(Labels2{j,k},'good') == 1
                 TimeStamps(j,k) = TimeStamps2(j,k);
             end
          end
        end
        
        indirect_unit = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        
        % Restructure to match input required
        indirect_cell = indirect_unit(~cellfun('isempty',indirect_unit));
        indirect_unit = cell2mat(indirect_unit(:));
        
        [Num_trials,T] = size(indirect_cell{1});
        units = length(indirect_cell);
        data_ = reshape(indirect_unit,Num_trials,units,T);
        
        
        for j=1:Num_trials
            D(j).data = squeeze(data_(j,:,:));
        end
        
        DataHigh(D,'DimReduce');
    end    
end       

%% Data high GPFA analysis  - M1

%%

clc;clear;



sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};

stroke= [4:5,9:10];

intact = [1:10];
intact = setdiff(intact, stroke);

Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';

tstart = {[0,0, 1, 6, 14, 5], [ 0, 0, 1, 6,  13, 5],...
    [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1],...
    [ 0,0,34, 8,52, 46,21,103,20,26],...
    [ 0,0,57, 7, 8, 50,42,55,9,20]...
    [ 0,0,0,12, 9, 5, 11,14]...
    [ 0,0,12, 60, 35, 44,2,9,19]...
    [0,0,0,0,1,23,1,24,1,21]...
    [0,0, 0,9,56,10 86, 20,1]...
    [0,0, 20,58,60,1]};

tstop  = {[0,0, 94, 106, 114, 106],[ 0, 0, 104, 105, 114, 105],...
    [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93],...
    [ 0,0,133, 68, 125, 86,103,151,262,109],...
    [ 0,0,128,123,72, 154,139,111,60,50]...
    [0,0, 0,200,120,105, 134,89]...
    [0,0, 74,83,84, 80,81,103,85]...
    [0,0,0,0,116,45,43,70,32,82]...
    [0,0, 0,109,88,82,145,101,150]...
    [0,0, 80,91,103,107]};

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6]};

% bmiBlocks = dir(rootpath); 
before_zero = 2.2;
after_zero =0.2;

tp = [];
tn = [];

res = [];
indirect = [];

% bin = [5,10,20,25,50,100]; % binnning at different ms
% for b=4%1:length(bin)
for s=1%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=4%robust_session{s}
    
        
    % Load timestamps
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
        bin_window = 1/1000; %in sec
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
        events = rewards_onset;
        valid_perf = performance(tstart{s}(n):tstop{s}(n)); 
                
        TimeStamps = cell(size(Labels1,1),size(Labels1,2));
        
        % get good trials only
        for j=1:size(Labels1,1)
          for k=2:size(Labels1,2)
             if strcmp(Labels1{j,k},'good') == 1
                 TimeStamps(j,k) = TimeStamps1(j,k);
             end
          end
        end
        
        indirect_unit = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        
        % Restructure to match input required
        indirect_cell = indirect_unit(~cellfun('isempty',indirect_unit));
        indirect_unit = cell2mat(indirect_unit(:));
        
        [Num_trials,T] = size(indirect_cell{1});
        units = length(indirect_cell);
        data_ = reshape(indirect_unit,Num_trials,units,T);
        
        
        for j=1:Num_trials
            D(j).data = squeeze(data_(j,:,:));
        end
        
        DataHigh(D,'DimReduce');
    end    
end       