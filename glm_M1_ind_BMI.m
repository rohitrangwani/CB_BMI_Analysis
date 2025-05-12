% GLM model - for M1 BMI
% Predictors were binned firing rates of either all Cb indirect units
% Response variable - M1 direct activity
% - Rohit Rangwani, 2023


%% Generate binned data

clc;clear;

Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Data_for_SpikeGLM\';

% bmiBlocks = dir(rootpath); 
% before_zero = 2.2;
% after_zero =0.2;

before_zero = 2.2;
after_zero =0.2;

tp = [];
tn = [];

res = [];
indirect = [];

bin = [5,10,25,50,100]; % binnning at different ms
for b=1:4%length(bin)

    bmiBlocks = dir('Z:\Data_for_SpikeGLM\');
    for n=4%:23
    
        
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
%         [tp_base,~,~] = fn_getPSTH_bar(TimeStamps_tp,events,Fs,savepath,2,0,trial_start,trial_stop);
        [tn_unit,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,before_zero,after_zero,trial_start,trial_stop);
%         [tn_base,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,2,0,trial_start,trial_stop); 
        
        tp_count = nnz(~cellfun(@isempty,tp_unit));
        tn_count = nnz(~cellfun(@isempty,tn_unit));

          if exist([rootpath,  bmiBlocks(n).name,'\Timestamps_B1.mat'],'file')
              load([rootpath,  bmiBlocks(n).name,'\Timestamps_B1.mat']);
              TimeStamps2 = TimeStamps2_B1;
              TimeStamps1 = TimeStamps1_B1;
              clear TimeStamps1_B1 TimeStamps2_B1;
          elseif exist([rootpath, bmiBlocks(n).name,'\Timestamps_B2.mat'],'file')
              load([rootpath,  bmiBlocks(n).name,'\Timestamps_B2.mat']);
              TimeStamps2 = TimeStamps2_B2;
              TimeStamps1 = TimeStamps1_B2;
              clear TimeStamps1_B2 TimeStamps2_B2;
          elseif exist([rootpath,  bmiBlocks(n).name,'\Timestamps.mat'],'file')
              load([rootpath,  bmiBlocks(n).name,'\Timestamps.mat']);
              TimeStamps2 = TimeStamps2_B1;
              TimeStamps1 = TimeStamps1_B1;
              clear TimeStamps1_B1 TimeStamps2_B1;
          else
              load([rootpath, bmiBlocks(n).name,'\Timestamps_Cb.mat']);        
              TimeStamps2(strcmp(Labels2,'good')~=0)=[];
              load([rootpath, bmiBlocks(n).name,'\Timestamps_M1.mat']);        
              TimeStamps1(strcmp(Labels1,'good')~=0)=[];
          end 
        
        tStamp = TimeStamps1(~cellfun(@isempty,TimeStamps1));  
        msize = length(tStamp);       
        tStamp = tStamp(randperm(msize, tp_count+tn_count));
        TimeStamps_tp = (tStamp(1:tp_count));
        TimeStamps_tn = (tStamp(tp_count+1:end));
        TimeStamps_tp = reshape(TimeStamps_tp,1,[]);
        TimeStamps_tn = reshape(TimeStamps_tn,1,[]);
        tp_idx = 1:tp_count;
        tn_idx = 1:tn_count;
        tp_unit = fn_getbinnedData(TimeStamps_tp,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,1,tp_idx);   
        
        if ~isempty(TimeStamps_tn)
            tn_unit = fn_getbinnedData(TimeStamps_tn,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,1,tn_idx);
%             tn_unit = cell2mat(tn_unit);
        end

        if ~tn_count
            res{b,n} = tp_unit;
        elseif ~tp_count
            res{b,n} = tn_unit;
        else
            res{b,n} = tp_unit - tn_unit;
        end   

        
%       indirect = [indirect;fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,trial_start(n),trial_stop(n),bin_window)];
        indirect_unit = fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,0,valid_perf);
        indirect{b,n} = indirect_unit;
       
        clear all_trial* TimeStamps2 tp_unit tn_unit
   end   

end
%% Reshape the data

indirect = reshape(indirect,size(indirect,1),size(indirect,2)*size(indirect,3));
res = reshape(res,size(res,1),size(res,2)*size(res,3));

%% Train the GLM model and test
  
  all_R2 = zeros(length(bin),10,10);
  
  r2 = zeros(length(bin),10,10);
  prob = []; 
  all_model = [];
  
  
  
%iterate for all binned samples  
for b=1:4%length(bin)    
    
% iterate for 10 time for crossvalidation    
for iter=1:10
    display(iter);

    res_cell  =  res(b,:);
    res_cell = res_cell(~cellfun('isempty',res_cell));
    indirect_cell  =  indirect(b,:);
    indirect_cell = indirect_cell(~cellfun('isempty',indirect_cell));
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
          for t=floor(dur/2):size(res_data,2)-(floor(dur/2)+1)
              start = start + 1;
              stop = stop + 1;

    %           start = start + 200/bin(b);
    %           stop = stop + 200/bin(b);
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
            SST = norm(z-mean(z))^2;
            R2 = 1 - SSE/SST;

            all_R2(b,iter,i) = R2;
        end
    end
end

end

%% Save the model/result and plot

save([rootpath, sub{s}, '\GLM_model_result_new.mat'], 'all_model','r2','all_R2', '-v7.3');
display("saved");

