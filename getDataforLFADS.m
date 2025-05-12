%% Get binnned data for lfads-torch
% Rohit 2023

tstart = {[0,0, 1, 6, 14, 5], [ 0, 0, 1, 6,  13, 5],...
    [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1],...
    [ 0,0,34, 8,52, 46,21,103,20,26],...
    [ 0,0,57, 7, 8, 50,42,55,9,20]...
    [ 0,0,0,12, 9, 5, 11,14]...
    [ 0,0,12, 60, 35, 44,2,9,19]...
    [0,0,0,0,1,23,1,24,1,21]...
    [0,0, 0,9,56,10 86, 20,1]...
    [0,0, 20,58,60,1]...
    [0,0,0, 66, 12,39,2, 1, 29,14, 19, 45,17,8,  2, 61, 4, 50, 38,29, 5, 1, 15],...
    [0,0,0, 9, 6, 11, 4, 45,19, 11, 6, 1, 20,43, 17,23]};

tstop  = {[0,0, 94, 106, 114, 106],[ 0, 0, 104, 105, 114, 105],...
    [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93],...
    [ 0,0,133, 68, 125, 86,103,151,262,109],...
    [ 0,0,128,123,72, 154,139,111,60,50]...
    [0,0, 0,200,120,105, 134,89]...
    [0,0, 74,83,84, 80,81,103,85]...
    [0,0,0,0,116,45,43,70,32,82]...
    [0,0, 0,109,88,82,145,101,150]...
    [0,0, 80,91,103,107]...
    [0,0, 0,109,68,61,10,50,83,127,88,101,55,80,65,117,47,159,76,90,80,35,57 ],...
    [0,0,0,85,80,93,100,82,140,55,100,38,60,120,60,82]};



sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160'};

stroke= [4:5,9:12];

intact = [1:12];
intact = setdiff(intact, stroke);
Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6],[8,14,15,16,18,19],[7,9,11,14]};
kl

before_zero = 4;
after_zero =0.0;

tp = [];
tn = [];

res = [];
indirect = [];


for s=12%2:12%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
    % Define save path
        savepath = [rootpath, '\lfads\Cb\', bmiBlocks(n).name,'_Cb.mat'];
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');
        events = rewards_onset;
        
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);

        valid_perf = performance(tstart{s}(n):tstop{s}(n)); 

        bin_window = 20/1000; % in seconds
  
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
%         events = rewards_onset;
        
        % Get only good units
%         TimeStamps = cell(size(Labels2,1),size(Labels2,2));
%         for j=1:size(Labels2,1)
%           for k=2:size(Labels2,2)
%              if strcmp(Labels2{j,k},'good') == 1
%                  TimeStamps(j,k) = TimeStamps2(j,k);
%              end
%           end
%         end
        
%       indirect = [indirect;fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
        indirect_unit = fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        indirect = indirect_unit(~cellfun('isempty',indirect_unit));
        n_units = size(indirect,1);
        [n_tr,n_time] = size(indirect{1});
        indirect = cell2mat(indirect);
        indirect = reshape(indirect,n_units,n_time,n_tr);
        events = logical(valid_perf < 15);
        
        save(savepath, 'indirect','events');
        
%         PD = 0.20 ;  % percentage 80%
%         cv = cvpartition(size(indirect,3),'HoldOut',PD);
%         train = indirect(:,:,cv.training);
%         test = indirect(:,:,cv.test);
        
        % Write ds to HDF5
%         filename = [rootpath,'\HDF5\',sub{s},'_',bmiBlocks(n).name];
%         h5write(filename,'train_encod_data',train);
%         h5write(filename,'valid_encod_data',test);
%         h5write(filename,'train_recon_data',train);
%         h5write(filename,'valid_recon_data',test);


   end   

end


%%
%% Get binnned data for lfads-torch

tstart = {[0,0, 1, 6, 14, 5], [ 0, 0, 1, 6,  13, 5],...
    [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1],...
    [ 0,0,34, 8,52, 46,21,103,20,26],...
    [ 0,0,57, 7, 8, 50,42,55,9,20]...
    [ 0,0,0,12, 9, 5, 11,14]...
    [ 0,0,12, 60, 35, 44,2,9,19]...
    [0,0,0,0,1,23,1,24,1,21]...
    [0,0, 0,9,56,10 86, 20,1]...
    [0,0, 20,58,60,1]...
    [0,0,0, 66, 12,39,2, 1, 29,14, 19, 45,17,8,  2, 61, 4, 50, 38,29, 5, 1, 15],...
    [0,0,0, 9, 6, 11, 4, 45,19, 11, 6, 1, 20,43, 17,23]};

tstop  = {[0,0, 94, 106, 114, 106],[ 0, 0, 104, 105, 114, 105],...
    [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93],...
    [ 0,0,133, 68, 125, 86,103,151,262,109],...
    [ 0,0,128,123,72, 154,139,111,60,50]...
    [0,0, 0,200,120,105, 134,89]...
    [0,0, 74,83,84, 80,81,103,85]...
    [0,0,0,0,116,45,43,70,32,82]...
    [0,0, 0,109,88,82,145,101,150]...
    [0,0, 80,91,103,107]...
    [0,0, 0,109,68,61,10,50,83,127,88,101,55,80,65,117,47,159,76,90,80,35,57 ],...
    [0,0,0,85,80,93,100,82,140,55,100,38,60,120,60,82]};



sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160'};

stroke= [4:5,9:12];

intact = [1:12];
intact = setdiff(intact, stroke);
Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6],[8,14,15,16,18,19],[7,9,11,14]};


before_zero = 2.2;
after_zero =0.2;

tp = [];
tn = [];

res = [];
indirect = [];


for s=2:12%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
    % Define save path
        savepath = [rootpath, '\lfads\M1\', bmiBlocks(n).name,'_M1.mat'];
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');
        events = rewards_onset;
        
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);

        valid_perf = performance(tstart{s}(n):tstop{s}(n)); 

        bin_window = 50/1000; % in seconds
  
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
%         events = rewards_onset;
        
        % Get only good units
%         TimeStamps = cell(size(Labels2,1),size(Labels2,2));
%         for j=1:size(Labels2,1)
%           for k=2:size(Labels2,2)
%              if strcmp(Labels2{j,k},'good') == 1
%                  TimeStamps(j,k) = TimeStamps2(j,k);
%              end
%           end
%         end
        
%       indirect = [indirect;fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
        indirect_unit = fn_getbinnedData(TimeStamps1,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        indirect = indirect_unit(~cellfun('isempty',indirect_unit));
        n_units = size(indirect,1);
        [n_tr,n_time] = size(indirect{1});
        indirect = cell2mat(indirect);
        indirect = reshape(indirect,n_units,n_time,n_tr);
        events = logical(valid_perf < 15);
        save(savepath, 'indirect','events');
        
%         PD = 0.20 ;  % percentage 80%
%         cv = cvpartition(size(indirect,3),'HoldOut',PD);
%         train = indirect(:,:,cv.training);
%         test = indirect(:,:,cv.test);
        
        % Write ds to HDF5
%         filename = [rootpath,'\HDF5\',sub{s},'_',bmiBlocks(n).name];
%         h5write(filename,'train_encod_data',train);
%         h5write(filename,'valid_encod_data',test);
%         h5write(filename,'train_recon_data',train);
%         h5write(filename,'valid_recon_data',test);


   end   

end

