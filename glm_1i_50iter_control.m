% GLM model - GLM-1i control analysis
% Predictors were binned firing rates of M1 indirects
% response variable were randomly selected and CB direct matched CB indirect units
% 50 iteration for cross-validation
% - Rohit Rangwani, 2023


%% Generate binned data

clc;clear;

%Load all session info
bmi_session_info;

before_zero = 2.2;
after_zero =0.2;

tp = [];
tn = [];

res = [];
indirect = [];

bin = [5,10,25,50,100]; % binnning at different ms
TU = 0;

for z=15

for b=4%1:4%length(bin)
for s=z%M1
    s
    rootpath = 'Z:\Rohit\BMI_Data\';
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}
    
        
    % Define save path
        savepath = [rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];
        load([rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Timestamps_Direct.mat']);
        load([rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');

        events = rewards_onset;
        load([rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Performance_stats_early_late.mat']);


        bin_window = bin(b)/1000; % in seconds

        
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
%         events = rewards_onset;
        
        % Get only good units
        TimeStamps = cell(size(Labels2,1),size(Labels2,2));
        for j=1:size(Labels2,1)
          for k=2:size(Labels2,2)
             if strcmp(Labels2{j,k},'good') == 1
                 TimeStamps(j,k) = TimeStamps2(j,k);
             end
          end
        end
        
        [tp_unit,~,~] = fn_getPSTH_bar(TimeStamps_tp,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
%         [tp_base,~,~] = fn_getPSTH_bar(TimeStamps_tp,events,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));
        [tn_unit,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
%         [tn_base,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n)); 
        
        tp_count = nnz(~cellfun(@isempty,tp_unit));
        tn_count = nnz(~cellfun(@isempty,tn_unit));
        



            [tp_unit,~,~] = fn_getPSTH_bar(TimeStamps,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
            [tp_base,~,~] = fn_getPSTH_bar(TimeStamps,events,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));
            [tn_unit,~,~] = fn_getPSTH_bar(TimeStamps,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
            [tn_base,~,~] = fn_getPSTH_bar(TimeStamps,events,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));    

            tp_unit =  cell2mat(tp_unit(:));
            tp_base =  cell2mat(tp_base(:));
            tn_unit =  cell2mat(tn_unit(:));
            tn_base =  cell2mat(tn_base(:));

            tStamp = TimeStamps(~cellfun(@isempty,TimeStamps));  

            TimeStamps_Tp = tStamp;
            TimeStamps_Tn = tStamp;            
        % Iterate to find 50 random response variables using M1 indirects
        for iter=(1:1000)
            clear tp_unit tn_unit           


            if ~TU

                TimeStamps_Tp_Iter = TimeStamps_Tp(randperm(length(TimeStamps_Tp), tp_count));
                TimeStamps_Tp_Iter = reshape(TimeStamps_Tp_Iter,1,[]);
                tp_idx = 1:tp_count;
                tp_unit = fn_getbinnedData(TimeStamps_Tp_Iter,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,1,tp_idx);   

                if tn_count              
                    TimeStamps_Tn_Iter = TimeStamps_Tn(randperm(length(TimeStamps_Tn), tn_count));
                    TimeStamps_Tn_Iter = reshape(TimeStamps_Tn_Iter,1,[]);
                    tn_idx = 1:tn_count;
                    tn_unit = fn_getbinnedData(TimeStamps_Tn_Iter,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,1,tn_idx);
                end

                if ~tn_count
                    res{b,n,iter} = tp_unit;
                elseif ~tp_count
                    res{b,n,iter} = tn_unit;
                else
                    res{b,n,iter} = tp_unit - tn_unit;
                end   
            
            else
                if(length(TimeStamps_TU) > tp_count+tn_count)
                    tu_count = tp_count+tn_count;
                else
                    tu_count = length(TimeStamps_TU);
                end    
                TimeStamps_TU_Iter = TimeStamps_TU(randperm(length(TimeStamps_TU), tu_count));
                TimeStamps_TU_Iter = reshape(TimeStamps_TU_Iter,1,[]);
                tu_idx = 1:tu_count;
                res{b,n,iter} = fn_getbinnedData(TimeStamps_TU_Iter,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,1,tu_idx);  
            end

        end
        
                load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
%         events = rewards_onset;
        
        % Get only good units
        TimeStamps = cell(size(Labels1,1),size(Labels1,2));
        for j=1:size(Labels1,1)
          for k=2:size(Labels1,2)
             if strcmp(Labels1{j,k},'good') == 1
                 TimeStamps(j,k) = TimeStamps1(j,k);
             end
          end
        end

        
        indirect_unit = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
                  
        indirect{b,n} = indirect_unit;
        
        
        clear all_trial* TimeStamps2 rewards_onse*
%         end
   end   
end
end

% Train the GLM model and test


  all_R2 = zeros(length(bin),50,20);
  r2 = zeros(length(bin),50,20);


  fileName = '\GLM_M1_Cb_control.mat';
  m = matfile([rootpath, sub{s},fileName],'Writable',true);
  m.r2 = r2;
  m.all_R2 = all_R2;
%   m.all_model = all_model;
  
%iterate for all binned samples  
for b=4%1:4%length(bin)    
    
% iterate for 10 time for crossvalidation    

for iter=1:1000

    disp(iter);
    res_cell  =  res(b,robust_session{s},iter);
    indirect_cell  =  indirect(b,robust_session{s});

    for i=1:length(indirect_cell)
        if(~isempty(res_cell{i}))

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
          
          shuffle = randperm(num_trial);
          for t=floor(dur/2):size(res_data,2)-(floor(dur/2)+1)
              start = start + 1;
              stop = stop + 1;

              rearranged = indirect_data(:,(start):(stop));
              rearranged = reshape(rearranged,num_trial,(dur+1)*num_indirect);
              %only 2 sec duration data eith different time lags
              
              % shuffle trials for control
              rearranged = rearranged(shuffle,:);
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
          if(iter==1)
              all_model{b,i} = mdl;
          end
          %R-squared value for model
          m.r2(b,iter,i) = mdl.Rsquared.Ordinary;

            Xnew = pred((9*end/10)+1:end,:);
            ypred = predict(mdl,Xnew);
            
            z = var((9*end/10)+1:end,:);
            r = z-ypred;
            normr = norm(r);
            SSE = normr.^2; 
            SST = norm(z-mean(z))^2;
            R2 = 1 - SSE/SST;

            m.all_R2(b,iter,i) = R2;

        end
        end
    end

end

end


%
 save([rootpath, sub{s},fileName],'all_model','-append');
 
end
