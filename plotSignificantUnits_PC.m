%% Generate binned data and significant units for purkinje 
%  - Rohit

clc;clear;

%Load all session info
bmi_session_info;

before_zero = 3;
after_zero =1;

tp = [];
tn = [];

res = [];
indirect = [];

for s=stroke % intact or stroke
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    load([rootpath,'\purkinje_units_all.mat'])
    load([rootpath,'\purkinje_ch_all.mat'])
    for n=intersection{s}
    
        
    % Define save path
        savepath = [rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];

        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');
        events = rewards_onset;
%         events = all_trials;
          % Define save path

        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);

        valid_perf = performance(tstart{s}(n):tstop{s}(n)); 
        
%       load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'valid_p*');
%       if exist('valid_pref','var') == 1
%         valid_perf = valid_pref;
%       end

%       remove_idx = (valid_perf==15);
%         events = rewards_onset(tstart{s}(n):tstop{s}(n));
%         events = events(valid_perf<15);
%         clear valid_perf
%         bin_window = bin(b)/1000; % in seconds
%         tp = [tp;fn_getbinnedData(TimeStamps_tp,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
%         tn = [tn;fn_getbinnedData(TimeStamps_tn,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
        
%         tp_unit = fn_getbinnedData(TimeStamps_tp,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,1,valid_perf);
%         tn_unit = fn_getbinnedData(TimeStamps_tn,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,1,valid_perf);
        [tp_unit,~,~] = fn_getPSTH_bar(TimeStamps_tp,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
        [tn_unit,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
        [tp_base,~,~] = fn_getPSTH_bar(TimeStamps_tp,all_trials,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));
        [tn_base,~,~] = fn_getPSTH_bar(TimeStamps_tn,all_trials,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));        
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
%         events = rewards_onset;
        
        % Get only PC units
        ch = getfield(channels,bmiBlocks(n).name);
        unit = getfield(units,bmiBlocks(n).name); 
        % Get only good units
        TimeStamps = cell(size(Labels2,1),size(Labels2,2));
        for j=1:size(Labels2,1)
          for k=2:size(Labels2,2)
             if strcmp(Labels2{j,k},'good') == 1
                 % Purkinje putative
                 idx = find(ch==k+1);
                 if ~isempty(idx) && unit(idx(1)) ~= 0
                    TimeStamps(j,k) = TimeStamps2(j,k);
                 end
                 % Not Purkinje 
%                  idx = find(ch~=k+1);
%                  if ~isempty(k) && unit(k(1)) ~= 0
%                     TimeStamps(j,k) = TimeStamps2(j,k);
%                  end
                
             end
          end
        end

%       indirect = [indirect;fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
  % Define save path
        savepath = [rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Indirect_Units_Reward_Cb_AR\'];
        [indirect_unit,~,~] = fn_getPSTH_bar(TimeStamps,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
        [indirect_base,~,~] = fn_getPSTH_bar(TimeStamps,all_trials,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));
%         indirect_unit = fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        
%         if isempty(tn_unit)
%             res{b,n} = tp_unit;
%         else
%             res{b,n} = tp_unit - tn_unit;
%         end       


%         end       
        tp{s,n} = tp_unit;
        tn{s,n} = tn_unit;
        indirect{s,n} = indirect_unit;
        base_tp{s,n} = tp_base;
        base_tn{s,n} = tn_base;
        base_indirect{s,n} = indirect_base;
%         save([rootpath, sub{s},'\Data\', bmiBlocks(n).name, '\GLM_model_data.mat'], 'indirect_unit','tp_unit','tn_unit');
        
%         indirect_ = indirect_unit(~cellfun('isempty',indirect_unit));
%         A = zeros(size(indirect_unit{1}));
%         for i=1:length(indirect_unit)
%             A = A+indirect_unit{i};
%         end
%         indirect{b,s,n} = A./(length(indirect_unit));

%       clear valid_pref valid_perf all_trial*
        
   end   

   

end

%%

indirect_ = reshape(indirect,size(indirect,1),size(indirect,2)*size(indirect,3));
tp_ = reshape(tp,size(tp,1),size(tp,2)*size(tp,3));
tn_ = reshape(tn,size(tn,1),size(tn,2)*size(tn,3));

indirect_b = reshape(base_indirect,size(base_indirect,1),size(base_indirect,2)*size(base_indirect,3));
tp_b = reshape(base_tp,size(base_tp,1),size(base_tp,2)*size(base_tp,3));
tn_b = reshape(base_tn,size(base_tn,1),size(base_tn,2)*size(base_tn,3));

tp_cell = tp_(~cellfun('isempty',tp_));
tn_cell = tn_(~cellfun('isempty',tn_));
indirect_cell = indirect_(~cellfun('isempty',indirect_));

tp_cell_b = tp_b(~cellfun('isempty',tp_b));
tn_cell_b = tn_b(~cellfun('isempty',tn_b));
indirect_cell_b = indirect_b(~cellfun('isempty',indirect_b));

% tp_cell = cell2mat(tp_cell);
% tn_cell = cell2mat(tn_cell);
tp_data = [];
tn_data = [];
indirect_data = [];

for i=1:length(indirect_cell)
    tp_data{i} = cell2mat(tp_cell{i}(:));
    tn_data{i} = cell2mat(tn_cell{i}(:));
    indirect_data{i} = cell2mat(indirect_cell{i}(:));
    b_tp_data{i} = cell2mat(tp_cell_b{i}(:));
    b_tn_data{i} = cell2mat(tn_cell_b{i}(:));
    b_indirect_data{i} = cell2mat(indirect_cell_b{i}(:));
end
indirect_data = cell2mat(indirect_data');

tp_data = cell2mat(tp_data');
tn_data = cell2mat(tn_data');

b_indirect_data = cell2mat(b_indirect_data');

b_tp_data = cell2mat(b_tp_data');
b_tn_data = cell2mat(b_tn_data');

%% Get significantly modulated units only ( 4 STD)

C = 4;
T = 800:3200;

idx1 = (C*std(b_indirect_data,[],2) + mean(b_indirect_data,2) < max(indirect_data(:,T),[],2));
idx2 = (mean(b_indirect_data,2) - C*std(b_indirect_data,[],2) > min(indirect_data(:,T),[],2));

% idx = idx1|idx2;
% 
% indirect_data = indirect_data(idx,:);
% indirect_data_M1(2.5*std(indirect_data_M1(:,1:(1/2)*end),[],2) + mean(indirect_data_M1(:,1:0.5*end),2) < max(indirect_data_M1(:,.5*end:end),[],2));

idx_tp = C*std(b_tp_data,[],2) + mean(b_tp_data,2) < max(tp_data(:,T),[],2);

if ~isempty(tn_data)
    idx_tn = ( mean(b_tn_data,2) - C*std(b_tn_data,[],2) > min(tn_data(:,T),[],2));
end

% idx_tp = 2*std(b_tp_data,[],2) + mean(b_tp_data,2) < max(tp_data(:,:),[],2);
% 
% if ~isempty(tn_data)
%     idx_tn = ( mean(b_tn_data,2) - 2*std(b_tn_data,[],2) > min(tn_data(:,:),[],2));
% end
%% Tp

figure;
tp_data_1 = tp_data(idx_tp,500:3500);
weights = tp_data_1(any(tp_data_1,2),:);
%         weight = normalize(weights,2,'norm',inf);
weight = normalize(weights,2,"range",[-1 1]);
[m,idx] = max(weight,[],2);
[ii,idx] = sort(idx); 

weight = weight(idx,:);
% weight = sortrows(weight,max(weight'));
image(weight, 'CDataMapping', 'scaled')
colorbar;
set(gcf,'WindowState','maximized')
set(gca, 'TickDir', 'out')


%% Tn
            figure;
if ~isempty(tn_data)
%             tn_data_1 = tn_data(idx_tn,:);
%             weights = tn_data_1(any(tn_data_1,2),:);

    weights = tn_data(any(tn_data,2),500:3500);
%             weight = normalize(weights,2,'norm',inf);
    weight = normalize(weights,2,"range",[-1 1]);
    [m,idx] = min(weight(:,:),[],2);
    [ii,idx] = sort(idx,'ascend'); 

    weight = weight(idx,:);
%   weight = weight([1:3,5:end],:);
    % weight = sortrows(weight,max(weight'));
    image(weight, 'CDataMapping', 'scaled')
    colorbar;
    set(gcf,'WindowState','maximized')
    set(gca, 'TickDir', 'out')
%             saveas(gcf,[rootpath,'\','10_TnUnits.tiff']);
%             close;
end        
%% Positive modulated significant indirect units
figure;
indirect_data1 = indirect_data(idx1,500:3500);
weights = indirect_data1(any(indirect_data1,2),:);
%         weight = zscore(weights,0,2);
weight = normalize(weights,2,"range",[-1 1]);
[m,idx] = max(weight,[],2);
[ii,idx] = sort(idx); 

weight = weight(idx,:);
% weight = sortrows(weight,max(weight'));
image(weight, 'CDataMapping', 'scaled')
colorbar;
set(gcf,'WindowState','maximized')
set(gca, 'TickDir', 'out')
%         saveas(gcf,[rootpath,'\','10_IndirectUnits_Max.tiff']);
%         close;
 %% Negative modulated significant indirect units
figure;
indirect_data2 = indirect_data(idx2,500:3500);
weights = indirect_data2(any(indirect_data2,2),:);
%         weight = zscore(weights,0,2);
weight = normalize(weights,2,"range",[-1 1]);
[m,idx] = min(weight(:,:),[],2);
[ii,idx] = sort(idx,'ascend'); 

weight = weight(idx,:);
% weight = sortrows(weight,max(weight'));
image(weight, 'CDataMapping', 'scaled')
colorbar;
set(gcf,'WindowState','maximized')
set(gca, 'TickDir', 'out')
%         saveas(gcf,[rootpath,'\','10_IndirectUnits_Min.tiff']);
%         close;       
%% All significant indirect units
figure;
idx = idx1 | idx2;
indirect_data1 = indirect_data(idx,500:3500);
weights = indirect_data1(any(indirect_data1,2),:);
weight = normalize(weights,2,"range",[-1 1]);
%         weight = zscore(weights,0,2);
%         weight = normalize(weights,2);
[m,idx] = max(weight(:,:),[],2);
[ii,idx] = sort(idx); 

weight = weight(idx,:);
% weight = sortrows(weight,max(weight'));
image(weight, 'CDataMapping', 'scaled')
colorbar;
set(gcf,'WindowState','maximized')
set(gca, 'TickDir', 'out')
%         saveas(gcf,[rootpath,'\','10SD_IndirectUnits.tiff']);
%         close;
           