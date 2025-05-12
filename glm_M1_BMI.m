% GLM model - for M1 BMI
% Predictors were binned firing rates of either all Cb indirect units
% Response variable - M1 direct activity
% - Rohit Rangwani, 2023

%% Generate binned data

clc;clear;

Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Data_for_SpikeGLM\';

before_zero = 0.2;
after_zero =2.2;

tp = [];
tn = [];

res = [];
indirect = [];

first = 4;
last = 23;

bin = [5,10,25,50,100]; % binnning at different ms
for b=2:4%1:4%length(bin)

    bmiBlocks = dir('Z:\Data_for_SpikeGLM\');
    for n=first:last
        disp(n);
        
    % Define save path
        savepath = [rootpath, bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];
        load([rootpath, bmiBlocks(n).name,'\Timestamps_Direct.mat']);
        load([rootpath, bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');
        if exist('all_trials_B1','var')
            events = all_trials_B1;
        else
            events = all_trials;
        end

        load([rootpath, bmiBlocks(n).name,'\Performance_stats_early_late.mat']);

        bin_window = bin(b)/1000; % in seconds
       
        [tp_unit,~,~] = fn_getPSTH_bar(TimeStamps_tp,events,Fs,savepath,before_zero,after_zero,trial_start,trial_stop);
        [tp_base,~,~] = fn_getPSTH_bar(TimeStamps_tp,events,Fs,savepath,2,0,trial_start,trial_stop);
        [tn_unit,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,before_zero,after_zero,trial_start,trial_stop);
        [tn_base,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,2,0,trial_start,trial_stop);    
        
        tp_unit =  cell2mat(tp_unit(:));
        tp_base =  cell2mat(tp_base(:));
        tn_unit =  cell2mat(tn_unit(:));
        tn_base =  cell2mat(tn_base(:));
        
        tp_idx = (1.5*std(tp_base,[],2) + mean(tp_base,2) < max(tp_unit,[],2));
        tp_unit = fn_getbinnedData(TimeStamps_tp,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,1,tp_idx);   
        
        if ~isempty(tn_unit)
            tn_idx = ( mean(tn_base,2) - 1.5*std(tn_base,[],2) > min(tn_unit,[],2));
            tn_unit = fn_getbinnedData(TimeStamps_tn,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,1,tn_idx);
        end

        if isempty(tn_unit)
            res{b,n} = tp_unit;
        elseif isempty(tp_unit)
            res{b,n} = tn_unit;
        else
            res{b,n} = tp_unit - tn_unit;
        end   

        if(~isempty(res{b,n}))
 
          if exist([rootpath,  bmiBlocks(n).name,'\Timestamps_B1.mat'],'file')
              load([rootpath,  bmiBlocks(n).name,'\Timestamps_B1.mat']);
              TimeStamps2 = TimeStamps2_B1;
              clear TimeStamps1_B1 TimeStamps2_B1;
          elseif exist([rootpath, bmiBlocks(n).name,'\Timestamps_B2.mat'],'file')
              load([rootpath,  bmiBlocks(n).name,'\Timestamps_B2.mat']);
              TimeStamps2 = TimeStamps2_B2;
              clear TimeStamps1_B2 TimeStamps2_B2;
          elseif exist([rootpath,  bmiBlocks(n).name,'\Timestamps.mat'],'file')
              load([rootpath,  bmiBlocks(n).name,'\Timestamps.mat']);
              TimeStamps2 = TimeStamps2_B1;
              clear TimeStamps1_B1 TimeStamps2_B1;
          else
              load([rootpath, bmiBlocks(n).name,'\Timestamps_Cb.mat']);        
%               TimeStamps2(strcmp(Labels2,'good')==0)=[];
              TimeStamps2(strcmp(Labels2,'good')~=0)=[];
          end 
          
        [id_unit,~,~] = fn_getPSTH_bar(TimeStamps2,events,Fs,savepath,before_zero,after_zero,trial_start,trial_stop);
        [id_base,~,~] = fn_getPSTH_bar(TimeStamps2,events,Fs,savepath,2,0,trial_start,trial_stop);  
       
        mean_id = cellfun(@mean,id_base);
        std_id = cellfun(@std,id_base);
%         min_id = cellfun(@min,id_unit,'UniformOutput',false);
        max_id = cellfun(@max,id_unit,'UniformOutput',false);
%         min_id(cellfun('isempty',min_id)) = {0};
        max_id(cellfun('isempty',max_id)) = {0};
%         id_idx1 = (2.5*std(id_base,[],2) + mean(id_base,2) < max(id_unit,[],2));
%         id_idx2 = ( mean(id_base,2) - 2.5*std(id_base,[],2) > min(id_base,[],2));
        max_id = cell2mat(max_id);
        id_idx1 = (3*std_id + mean_id > max_id);
%         id_idx2 = ( mean_d - 2.5*std_id > min_id);
%         idx = id_idx1 | id_idx2;
        idx = id_idx1;

        
%       indirect = [indirect;fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,trial_start(n),trial_stop(n),bin_window)];
        indirect_unit = fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,0,idx);         
        indirect{b,n} = indirect_unit;

        clear all_trial* TimeStamps2 rewards_onset*
        else
            disp("No significant response variable for n:");
            disp(n);
        end
        
   end   

end

%% Train the GLM model and test

  all_model = [];

%iterate for all binned samples  
for b=4%2:4%1:4%length(bin)    
    
% iterate for 10 time for crossvalidation    
for iter=1:10
    display(iter);

    res_cell  =  res(b,first:last);
    indirect_cell  =  indirect(b,first:last);
    indirect_data = [];

    res_data = [];
    
    for i=1:length(indirect_cell)
        indirect_data = [];

        res_data = [];

        indirect_data1 = cell2mat(indirect_cell{i}(:));
        res_data = cell2mat(res_cell(i));
        num_trial = size(res_data, 1); 

        if(~isempty(indirect_data1))
            
        indirect_data =  indirect_data1;
        num_indirect =floor(size(indirect_data, 1)/size(res_data, 1));

          start = 0;
          
%           dur = 400/bin(b);
          dur = 100/bin(b);
%           stop = dur -1; 
          stop = dur; 
          pred = [];
          var = [];
          rearranged = [];
          for t=floor(dur/2)+1:size(res_data,2)-(floor(dur/2)+1)
              start = start + 1;
              stop = stop + 1;

              rearranged = indirect_data(:,(start):(stop));
              rearranged = reshape(rearranged,num_trial,(dur+1)*num_indirect);
              %only 2 sec duration data eith different time lags
              pred{t} = rearranged;

    %           var = resp(:,start:stop); %response variable
              var{t} = res_data(:,t); %% response variable is the activity at reward(2s)

          end
%           mdl = fitglm(pred(floor(1:9*end/10),:), var);
          pred = cell2mat(pred');
          var = cell2mat(var');
          rand_idx = randperm(size(var, 1));
          var = var(rand_idx,:);
          pred = pred(rand_idx,:);
          mdl = fitglm(pred(1:(9*end/10),:), var(1:(9*end/10)));  %% response variable is the data at reward
%           mdl = fitglm(pred(1:(8*end/10),:), var(1:(8*end/10)));  %% response variable is the data at reward
          all_model{b,iter,i} = mdl;
          
          %R-squared value for model
          r2(b,iter,i) = mdl.Rsquared.Ordinary;
          adj_r2(b,iter,i) =  mdl.Rsquared.Adjusted;
          prob{b,iter,i} = mdl.Fitted.Response;
          
            Xnew = pred((9*end/10)+1:end,:);
            ypred = predict(mdl,Xnew);
            z = var((9*end/10)+1:end,:);
            r = z-ypred;
            normr = norm(r);
            SSE = normr.^2; 
            SST = norm(z-mean(z))^2;;
            R2 = 1 - SSE/SST;

            all_R2(b,iter,i) = R2;
        end
    end
end

end

%% Save the model/result and plot

save([rootpath, '\GLM_model_result_M1d_CbId_TaskStart_allData.mat'], 'all_model','r2','all_R2', '-v7.3');
display("saved"); 
