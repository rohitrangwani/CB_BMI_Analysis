% GLM model - for M1 predicting CB direct - GLM1-d
% Predictors were binned firing rates of M1 indirects
% Rohit 2023


%% Generate binned data

% clc;
clear;

%Load all session info
bmi_session_info;

before_zero = 2.2;
after_zero =0.2;

tp = [];
tn = [];

res = [];
indirect = [];

bin = [5,10,25,50,100]; % binnning at different ms

% loop for all the bins and all the subjects
for b=2:4%length(bin)
for s=15%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));

    %loop over all the robust sessions
    for n=robust_session{s}
    
        
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
%         
%         load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1_new.mat']);

        
        % Get only good units
        TimeStamps = cell(size(Labels1,1),size(Labels1,2));
        for j=1:size(Labels1,1)
          for k=2:size(Labels1,2)
             if strcmp(Labels1{j,k},'good') == 1
                 TimeStamps(j,k) = TimeStamps1(j,k);
             end
          end
        end
        
        % Get indirect binned data
        indirect_unit = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        
        if isempty(tn_unit)
            res{b,n} = tp_unit;
        elseif isempty(tp_unit)
            res{b,n} = tn_unit;
        else
            res{b,n} = tp_unit - tn_unit;
        end       
        
        indirect{b,n} = indirect_unit;
        

        end
   end   

end
end
%% Reshape the data

% indirect = reshape(indirect,size(indirect,1),size(indirect,2)*size(indirect,3));
% 
% res = reshape(res,size(res,1),size(res,2)*size(res,3));

%% Train the GLM model and test
   
  all_R2 = zeros(length(bin),10,10);
  
  r2 = zeros(length(bin),10,10);
  prob = []; 
  all_model = [];
  
  
  
%iterate for all binned samples  
for b=2:4%2:length(bin)    

% iterate for 10 time for crossvalidation    

for iter=1:10
    display(iter);


    res_cell  =  res(b,:);
    res_cell = res_cell(~cellfun('isempty',res_cell));

    % organize indirect data
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
      dur = 100/bin(b);
      stop = dur; 
      pred = [];
      var = [];
      for t=floor(dur/2):size(res_data,2)-(floor(dur/2)+1)
          start = start + 1;
          stop = stop + 1;

          rearranged = indirect_data(:,(start):(stop));
          rearranged = reshape(rearranged,num_trial,(dur+1)*num_indirect);
          %only 2 sec duration data eith different time lags
          pred{t} = rearranged;

          var{t} = res_data(:,t); %% response variable is the activity at reward(2s)

      end

          pred = cell2mat(pred');
          var = cell2mat(var');
          rand_idx = randperm(size(var, 1));
          var = var(rand_idx,:);
          pred = pred(rand_idx,:);
          mdl = fitglm(pred(1:(9*end/10),:), var(1:(9*end/10)));  %% response variable is the data at reward

          all_model{b,iter,i} = mdl;
          
          %R-squared value from model
          r2(b,iter,i) = mdl.Rsquared.Ordinary;
          adj_r2(b,iter,i) =  mdl.Rsquared.Adjusted;
          prob{b,iter,i} = mdl.Fitted.Response;
          
            Xnew = pred((9*end/10)+1:end,:);
            ypred = predict(mdl,Xnew);
            
            % calculating R2 for verfication
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

%% Save the model/result and plot

save([rootpath, sub{s}, '\GLM_model_M1_result_new.mat'], 'all_model','r2','all_R2', '-v7.3');


%% Get indirect weights and probability for all sessions

bmi_session_info;

for s=M1
%     clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat']);
    for j=1:size(all_model,3)
    all_bins = [];
    for k=1:10
        figure;
        model = all_model{2,k,j};
        weights = model.Coefficients.Estimate(2:end);
        weights = reshape(abs(weights),length(weights)/11,11);

        weights = weights(any(weights,2),:);


        weight = normalize(weights,2,'norm',inf);

        [m,idx] = nanmax(weight,[],2);
        [ii,idx] = sort(idx); 
        

        weight = weight(idx,:);
        % weight = sortrows(weight,max(weight'));
        image(weight, 'CDataMapping', 'scaled')
        colorbar;
%         set(gcf,'WindowState','maximized')

        saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_M1_Sesssion_new_',num2str(j),'_',num2str(k),'.fig']);
        close;
        all_bins{k}=ii;
    end
    figure;
    [GC,GR] = groupcounts(cell2mat(all_bins(:))); 
%     histogram(cell2mat(all_bins(:)));
    bar(GC/10.0);
    saveas(gcf,[rootpath,'\',sub{s},'\','Histogram_M1_Sesssion__new_',num2str(j),'.fig']);
    close;
    end


    for j=1:size(all_R2,3)
        avg = mean((all_R2),2);
        subplot(5,2,j);
        bar(avg(:,:,j));
    end
    saveas(gcf,[rootpath,'\',sub{s},'\','M1_Test_R2_new_.tiff']);
    close;
    figure;
    for j=1:size(all_R2,3)
        avg = mean((r2),2);
        subplot(5,2,j);
        bar(avg(:,:,j));
    end
    saveas(gcf,[rootpath,'\',sub{s},'\','M1_Training_R2_new_.tiff']);
    close;
end


%% Get highest weight time lag histogram

bmi_session_info;

bin = [5,10,25,50,100];

str_bins = [];
int_bins = [];

for s=M1_stroke
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat']);
    for j=1:size(all_model,3) 
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
        for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,k,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

%                 set(gcf,'WindowState','maximized')
% 
%                 saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                str_bins{s,j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
        end
        end
    end
end 

for s=M1_intact
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat']);
    for j=1:size(all_model,3)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
        for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,k,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

%                 set(gcf,'WindowState','maximized')
% 
%                 saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                int_bins{s,j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
        end
        end
    end
end
%% plot histogram with kernel density estimates - intact
    figure;
    
    bin_2 = int_bins(:,:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    yline(mean(counts)); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;
     for i=2:4

        temp = int_bins(:,:,:,i);
        tmp = max(cell2mat(temp(:)));
        pd = fitdist(cell2mat(temp(:)).*(12/(tmp+1)),'Kernel');

        xlim([0,12]);
        x = 0.5:.05:12;
        ySix = pdf(pd,x);
        plot(x,ySix); hold on;
    end
    legend('uniform', 'histogram', '10 ms', '25 ms', '50 ms');  
%% plot histogram with kernel density estimates - stroke
    figure;
    
    bin_2 = str_bins(:,:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    yline(mean(counts)); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;
     for i=2:4

        temp = str_bins(:,:,:,i);
        tmp = max(cell2mat(temp(:)));
        pd = fitdist(cell2mat(temp(:)).*(12/(tmp+1)),'Kernel');

        xlim([0,12]);
        x = 0.5:.05:12;
        ySix = pdf(pd,x);
        plot(x,ySix); hold on;
    end
    legend('uniform', 'histogram', '10 ms', '25 ms', '50 ms');   
%%

    figure;
    
    [GC,GR] = groupcounts(cell2mat(all_bins(:))); 
%     histogram(cell2mat(all_bins(:)));
    bar(GC/10.0);
%     saveas(gcf,[rootpath,'\','Histogram_M1_new.tiff']);
%     close;
%%    
 figure
 pd = fitdist(cell2mat(temp(:)),'Kernel','Kernel','epanechnikov');

xlim([0,4]);
x = 0:.05:4;
ySix = pdf(pd,x);
plot(x,ySix,'k-','LineWidth',2)
%%
    figure;
    
    counts = histcounts(cell2mat(all_bins(:,:,:,2)),'Normalization','probability');
    yline(mean(counts)); hold on;
    bar(counts);
%%
    figure;
    
    bin_2 = all_bins(:,:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
%     yline(mean(counts)); hold on;
    ylim([0,.25]);
    bar(counts); hold on;
    
 %% ks test
 str_bin = str_bins(:,:,:,2);
 str_counts = histcounts(cell2mat(str_bin(:)),'Normalization','probability');
 int_bin = int_bins(:,:,:,2);
 int_counts = histcounts(cell2mat(int_bin(:)),'Normalization','probability');
 figure;
%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(str_bin(:)), cell2mat(int_bin(:)), 10000, 'plotresult',1)
 %%
 bar(int_counts); hold on;
 bar(str_counts); hold on;
%% check for uniformity

[h,p,ksstat] = kstest2(int_counts,repelem(mean(int_counts),11))
[h,p,ksstat] = kstest2(str_counts,repelem(mean(str_counts),11))
[h,p,ksstat] = kstest2(str_counts,int_counts)

%%
uniform = repelem([1:11],(length(cell2mat(str_bin(:)))/11));
[h,p,ksstat] = kstest2(cell2mat(int_bin(:))',uniform)
[h,p,ksstat] = kstest2(cell2mat(str_bin(:))',uniform)
[h,p,ksstat] = kstest2(cell2mat(str_bin(:)), cell2mat(int_bin(:)))
%% bootstrap
figure;
[ci,bootstat] = bootci(10000,@kstest2,int_counts,repelem(mean(int_counts),11));
ci
hist(bootstat); hold on;
vline(ci(2)); hold on;
vline(immse(int_counts,repelem(mean(int_counts),11)),'k');
[h,p]= ttest(bootstat)
figure;
[ci,bootstat] = bootci(10000,@kstest2,str_counts,repelem(mean(str_counts),11));
ci
hist(bootstat); hold on;
vline(ci(2)); hold on;
vline(immse(str_counts,repelem(mean(str_counts),11)),'k');
[h,p]= ttest(bootstat)

%%
figure;
[ci,bootstat] = bootci(10000,@kstest2,str_counts,int_counts);
ci
hist(bootstat); hold on;
vline(ci(1)); hold on;
vline(ci(2)); hold on;
vline(immse(str_counts,int_counts),'k');
[h,p]= ttest(bootstat)

%%
samples = datasample(str_counts,10000);
samples2 = datasample(int_counts,10000);

[h,p,ksstat] = kstest2(samples,samples2)


%% R2 healthy vs stroke - for M1 
bmi_session_info;

% 4 sessions are not significant glm sessions
% sig_session = {[4,5],[3,4,6],[7,8,9,10,13,14],[3,4,5,8,9,10],[5,7,9],...
%     [4,6,8],[5:9],[8],[8],[3,6],[8,14,15,16,18]};
% sig_session = {[1,1,0,0,0,0,0,0,0,0],[1,1,1,0,0,0,0,0,0,0],[0,1,1,1,1,1,1,0,0,0],[1,1,1,1,1,1,0,0,0,0],[1,0,1,1,0,0,0,0,0,0],...
%     [1,1,1,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[1,0,1,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0]};

val_healthy = [];
val_stroke = [];
for s=M1_intact
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat'],'all_R2');
    
%     load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat'],'all_R2','all_model');
%     healthy_mdl = all_model;
    for j=1:size(all_R2,3)
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        
        if(r2(4))
            val_healthy{s,j} = r2(4);    
        end
    end
    end
end 

for s=M1_stroke
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat'],'all_R2');
    
%     load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat'],'all_R2','all_model');
%     stroke_mdl = all_model;
    for j=1:size(all_R2,3)
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(4))
            val_stroke{s,j} = r2(4);    
        end
    end
    end
end 

val1 = cell2mat(val_healthy(:));

val2 = cell2mat(val_stroke(:));
%%
figure('Color','white'); hold all;
bar([mean(val1),mean(val2)],0.6,'w');
errorbar([1,2],[mean(val1),mean(val2)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1))]);
scatter([repelem(1,length(val1)),repelem(2,length(val2))],[val1;val2],'k');
[h,p] = ttest2(val1,val2)
set(gca,'TickDir','out');
xlim([0.5,2.5]);
ylim([0,0.7]);
%% Violinplot
figure;
val2 = [val2; nan; nan]; % append NaN at end to use violin plot
violinplot([val1,val2]); hold on;
% plot([val1';val2'],'k');
set(gca,'TickDir','out');
%% Mixed effect analysis

rat_id = []
for col=1:(size(val_healthy,2))
    for row=1:(size(val_healthy,1))
    if(val_healthy{row,col})
        rat_id = [rat_id,row];
    end
    end
end
rat_id
for col=1:(size(val_stroke,2))
    for row=1:(size(val_stroke,1))
    if(val_stroke{row,col})
        rat_id = [rat_id,row+size(val_healthy,2)];
    end
    end
end
rat_id

tbl = [[val1(:,1);val2(:,1)],[zeros(length(val1),1);...
    ones(length(val2),1)],rat_id'];
tbl = array2table(tbl,'VariableNames',{'R2','IS','RatID'});
formula = 'R2 ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)


%% R2 Cb vs M1 

bmi_session_info;

val_healthy = [];
val_M1 = [];
val_M1_Cb = [];
clear val1 val2 val3;
for s=M1_stroke %M1_intact
    clear all_model all_R2 r2_Cb ;
    load([rootpath, sub{s}, '\GLM_model_result_new.mat'],'all_R2');
    for j=1:size(all_R2,3)
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2_Cb = avg(:,:,j);
        if(r2_Cb(4))
            val_healthy{s,j} = r2_Cb(4);    
        end
    end
    end
    clear all_model all_R2 r2_M1;
    load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat'],'all_R2');
    for j=1:size(all_R2,3)
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2_M1 = avg(:,:,j);
        if(r2_M1(4))
            val_M1{s,j} = r2_M1(4);    
        end
    end
    end
    clear all_model all_R2 r2_M1;
    load([rootpath, sub{s}, '\GLM_M1_Cb.mat'],'all_R2');
    for j=1:10
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2_M1 = avg(:,:,j);
        if(r2_M1(4))
            val_M1_Cb{s,j} = r2_M1(4);    
        end
    end
    end
end 

val1 = cell2mat(val_healthy(:));
val2 = cell2mat(val_M1(:));
val3 = cell2mat(val_M1_Cb(:));
p_val = ranksum(val1,val2,'tail','right');

[h,p] = ttest(val1,val2,'tail','right');

%% compare GLM-id vs GLM-1d
figure('Color','white'); hold all;
bar([mean(val1),mean(val2)],0.6,'w');
errorbar([1,2],[mean(val1),mean(val2)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1))]);
plot([val1';val2'],'k');
set(gca,'TickDir','out');
xlim([0.5,2.5]);
ylim([0,0.7]);
%     saveas(gcf,[rootpath,'\','R2-CbvsM1.tiff']);
%     close;
%% Mixed effect analysis, GLM-id vs GLM-1d

rat_id = []
for col=1:(size(val_healthy,2))
    for row=1:(size(val_healthy,1))
    if(val_healthy{row,col})
        rat_id = [rat_id,row];
    end
    end
end
rat_id
for col=1:(size(val_M1,2))
    for row=1:(size(val_M1,1))
    if(val_M1{row,col})
        rat_id = [rat_id,row+size(val_healthy,2)];
    end
    end
end
rat_id

tbl = [[val1(:,1);val2(:,1)],[zeros(length(val1),1);...
    ones(length(val2),1)],rat_id'];
tbl = array2table(tbl,'VariableNames',{'R2','IS','RatID'});
formula = 'R2 ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)
%%
figure;
violinplot([val1,val2]);
plot([val1';val2'],'k');
set(gca,'TickDir','out');
%% compare GLM-1i vs GLM-1d
figure('Color','white'); hold all;
bar([mean(val3),mean(val2)],0.6,'w');
errorbar([1,2],[mean(val3),mean(val2)],[std(val3)/sqrt((size(val3,1)-1)),std(val2)/sqrt((size(val2,1)-1))]);
plot([val3';val2'],'k');

[h,p] = ttest(val3,val2,'tail','right')
set(gca,'TickDir','out');
xlim([0.5,2.5]);
ylim([0,0.7]);

%% compare GLM-id vs GLM-1d vs GLM-1i
figure('Color','white'); hold all;
bar([mean(val1),mean(val2),mean(val3)],0.6,'w');
errorbar([1,2,3],[mean(val1),mean(val2),mean(val3)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1)),std(val3)/sqrt((size(val3,1)-1))]);
plot([val1';val2';val3'],'k');
set(gca,'TickDir','out');
xlim([0.5,3.5]);
ylim([0,0.7]);


%% Mixed effect analysis, compare GLM-1i vs GLM-1d

rat_id = []
for col=1:(size(val_M1_Cb,2))
    for row=1:(size(val_M1_Cb,1))
    if(val_M1_Cb{row,col})
        rat_id = [rat_id,row];
    end
    end
end
rat_id
for col=1:(size(val_M1,2))
    for row=1:(size(val_M1,1))
    if(val_M1{row,col})
        rat_id = [rat_id,row+size(val_M1_Cb,2)];
    end
    end
end
rat_id

tbl = [[val3(:,1);val2(:,1)],[zeros(length(val3),1);...
    ones(length(val2),1)],rat_id'];
tbl = array2table(tbl,'VariableNames',{'R2','IS','RatID'});
formula = 'R2 ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)
%%
figure;
v = [1 2];
groups = repelem(v,[length(val1)  length(val2)]);
val = [val1;val2];
violinplot(val,groups);

% ttest
[h,p_val] = ttest2(val1,val2,'Tail','right','Vartype','unequal');


%% PLot all R2 for all session

bmi_session_info;

val_healthy = [];
val_stroke = [];
figure;
xlim([0,4]); hold on;
xticks([0,1,2,3,4]);
xticklabels({'','10 ms','25 ms', '50 ms',''});
ylabel('R2');
count = 1;
r2 =[];

for s=M1_intact
    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_M1_result_new.mat'],'all_R2');
    for j=1:size(all_R2,3)
        if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2{count} = avg(:,:,j);
        tmp = cell2mat(r2(count));
        if(sum(tmp))
            plot(tmp(2:4),'.-','MarkerSize',15,'LineWidth',2); hold on;
            count = count +1;
        end
        end
    end

end

yline(0);
%% Percentage of unsucessful change vs R2
% clc;clear;
% rootpath = 'Z:\Rohit\BMI_Data\I112\Data\';
rootpath = 'Z:\Rohit\BMI_Data\';
sub = {'I096','I107','I110','I111','I112'};
first = [4,3,4,4,3];
last = [6,6,10,10,10];
bmiBlocks = dir(rootpath);
for s=1:5
    bmiBlocks = dir([rootpath, sub{s},'\Data\']);  
    for i=first(s):last(s)

        load([rootpath, sub{s},'\Data\',bmiBlocks(i).name,'\Performance_stats_early_late.mat'],'valid_p*');

        load([rootpath, '\', sub{s} '\GLM_model_result.mat'],'all_R2');
        avg = mean((all_R2),2);
        display(i-2);
        r2 = avg(:,:,i-(first(s)-1));
        %   subplot(5,2,j);
        if exist('valid_pref','var') == 1
        valid_perf = valid_pref;
        end

        early = valid_perf(1:floor(length(valid_perf)/3));
        late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 

        tOutcome(i,1) = sum(early==15)/sum(length(early))*100 - sum(late==15)/sum(length(late))*100;
        % successful late  % successful early
        tOutcome(i,2) = r2(4);
        clear valid_pref valid_perf r2


    end
    c = linspace(1,10,length(tOutcome));
    scatter(tOutcome(:,1),tOutcome(:,2),40,c,'filled')
    ax = gca;
    ax.XAxisLocation = 'origin'
    ax.YAxisLocation = 'origin'

    saveas(gcf,[rootpath,'\',sub{s},'\','Perf_vs_R2.tiff']);
    clear tOutcome
    close;
end
%%
figure('Color','white'); hold all;
bar(mean(tOutcome));
errorbar([1,2],mean(tOutcome),std(tOutcome)/sqrt((size(tOutcome,1)-1)));
% scatter(ones(size(tOutcome,1),1),tOutcome(:,1));
% scatter(ones(size(tOutcome,1),1)*2,tOutcome(:,2));
plot(tOutcome','k');
ylim([0 70]);
[~,pV] = ttest(tOutcome(:,1),tOutcome(:,2));
display(pV);

