% GLM model - GLM-id control analysis
% Predictors were binned firing rates of either all Cb indirect units
% Response variable - Cb direct activity
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
for b=4%2:4%length(bin)
for s=15%:10%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
    % Define save path
        savepath = [rootpath, sub{s},'\Data\',bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');
        events = rewards_onset;
        
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);

        valid_perf = performance(tstart{s}(n):tstop{s}(n)); 

        bin_window = bin(b)/1000; % in seconds
        
        [tp_unit,~,~] = fn_getPSTH_bar(TimeStamps_tp,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
        [tp_base,~,~] = fn_getPSTH_bar(TimeStamps_tp,all_trials,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));
        [tn_unit,~,~] = fn_getPSTH_bar(TimeStamps_tn,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
        [tn_base,~,~] = fn_getPSTH_bar(TimeStamps_tn,all_trials,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));    
        
        tp_unit =  cell2mat(tp_unit(:));
        tp_base =  cell2mat(tp_base(:));
        tn_unit =  cell2mat(tn_unit(:));
        tn_base =  cell2mat(tn_base(:));
        
        tp_idx = (2*std(tp_base,[],2) + mean(tp_base,2) < max(tp_unit,[],2));
        tp_unit = fn_getbinnedData(TimeStamps_tp,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,1,tp_idx);   
        
        if ~isempty(tn_unit)
            tn_idx = ( mean(tn_base,2) - 2*std(tn_base,[],2) > min(tn_base,[],2));
            tn_unit = fn_getbinnedData(TimeStamps_tn,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,1,tn_idx);
        end

        if isempty(tn_unit)
            res{b,n} = tp_unit;
        elseif isempty(tp_unit)
            res{b,n} = tn_unit;
        else
            res{b,n} = tp_unit - tn_unit;
        end   

        if(~isempty(res{b,n}))
            
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
        
%       indirect = [indirect;fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
        indirect_unit = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
                  
        indirect{b,n} = indirect_unit;

        end
   end   

end
end

% Train the GLM model and test
   
  all_R2 = zeros(length(bin),10,10);
  
  r2 = zeros(length(bin),10,10);
  prob = []; 
  all_model = [];

%iterate for all binned samples  
for b=4%2:4%length(bin)    
    
% iterate for 10 time for crossvalidation    
for iter=1:1000
%     display(iter);

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

        
        indirect_data =  indirect_data1;
        num_indirect =floor(size(indirect_data, 1)/size(res_data, 1));
      
          start = 0;
          
%           dur = 400/bin(b);
          dur = 100/bin(b);
%           stop = dur -1; 
          stop = dur; 
          pred = [];
          var = [];
          
          shuffle = randperm(num_trial);
          for t=floor(dur/2):size(res_data,2)-(floor(dur/2)+1)
              start = start + 1;
              stop = stop + 1;

              rearranged = indirect_data(:,(start):(stop));
              rearranged = reshape(rearranged,num_trial,(dur+1)*num_indirect);

              % shuffle trials for control
              rearranged = rearranged(shuffle,:);
              pred{t} = rearranged;
              var{t} = res_data(:,t); %% response variable is the activity at reward(2s)

          end
          pred = cell2mat(pred');
          var = cell2mat(var');
          rand_idx = randperm(size(var, 1));
          var = var(rand_idx,:);
          pred = pred(rand_idx,:);
          mdl = fitglm(pred(1:(9*end/10),:), var(1:(9*end/10)));  %% response variable is the data at reward
%           mdl = fitglm(pred(1:(8*end/10),:), var(1:(8*end/10)));  %% response variable is the data at reward
          
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

% Save the model/result and plot

save([rootpath, sub{s}, '\GLM_model_result_control_new.mat'], 'r2','all_R2', '-v7.3');
display("saved");



%% Compare R2 with control

rootpath = 'Z:\Rohit\BMI_Data\';
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160','I161','I170','I172'};
val_healthy = [];
val_stroke = [];
r2  = {};
figure;
xlim([0,4]); hold on;
xticks([0,1,2,3,4]);
xticklabels({'','10 ms','25 ms', '50 ms',''});
ylabel('R2');
count = 0;
tmp_1 = [];
tmp_2 = [];
for s=15%1:10

    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result.mat'],'all_R2');
    for j=1:size(all_R2,3)
        avg = mean((all_R2),2);
        r2{s,j} = avg(:,:,j);
        tmp = cell2mat(r2(s,j));
        tmp_1 = [tmp_1,tmp(4)];

    end
    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result_control_new.mat'],'all_R2');
    for j=1:size(all_R2,3)
        avg = mean((all_R2),2);
        r2{s,j} = avg(:,:,j);
        tmp = cell2mat(r2(s,j));
        tmp_2 = [tmp_2,tmp(4)];
        if(sum(tmp))
            plot(tmp(4),'.-','MarkerSize',15,'LineWidth',2); hold on;
            count = count +1;
        end
    end

end 
% yline(0);
% 
% disp(tmp_1(4) > tmp_2(4));

%% Compare R2 with control percentile

rootpath = 'Z:\Rohit\BMI_Data\';

val_healthy = [];
val_stroke = [];
r2  = {};

count = 0;
tmp_1 = [];
tmp_2 = [];
for s=15%:13

    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result.mat'],'all_R2');
    for j=1:size(all_R2,3)
        avg = mean((all_R2),2);
        r2{s,j} = avg(:,:,j);
        tmp = cell2mat(r2(s,j));
        tmp_1 = [tmp_1,tmp(4)];

    end
    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result_control_new.mat'],'all_R2');
    for j=1:size(all_R2,3)
        avg = mean((all_R2),2);
        r2{s,j} = avg(:,:,j);
        tmp = cell2mat(r2(s,j));
        tmp_2  = [tmp_2,prctile(all_R2(4,:,j),95)];
        
    end
end 
% yline(0);
% 
(tmp_1 > tmp_2)
