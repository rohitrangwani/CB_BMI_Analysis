% GLM model - for M1 BMI
% Predictors were binned firing rates of either all Cb indirect units
% Response variable - M1 indirect activity, 50 iter
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
TU = 1;

first = 19;
last =38;
bin = [5,10,25,50,100]; % binnning at different ms
for b=2:4%1:4%length(bin)

    bmiBlocks = dir('Z:\Data_for_SpikeGLM\');
    for n=first:last
    
        
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
        

            [tp_unit,~,~] = fn_getPSTH_bar(TimeStamps1,events,Fs,savepath,before_zero,after_zero,trial_start,trial_stop);
            [tp_base,~,~] = fn_getPSTH_bar(TimeStamps1,events,Fs,savepath,2,0,trial_start,trial_stop);
            [tn_unit,~,~] = fn_getPSTH_bar(TimeStamps1,events,Fs,savepath,before_zero,after_zero,trial_start,trial_stop);
            [tn_base,~,~] = fn_getPSTH_bar(TimeStamps1,events,Fs,savepath,2,0,trial_start,trial_stop);    

            tp_unit =  cell2mat(tp_unit(:));
            tp_base =  cell2mat(tp_base(:));
            tn_unit =  cell2mat(tn_unit(:));
            tn_base =  cell2mat(tn_base(:));

            tStamp = TimeStamps1(~cellfun(@isempty,TimeStamps1));  
            tp_idx = (3*std(tp_base,[],2) + mean(tp_base,2) < max(tp_unit,[],2));
            tn_idx = ( mean(tn_base,2) - 3*std(tn_base,[],2) > min(tn_unit,[],2));
            TimeStamps_Tp = tStamp(tp_idx);
            TimeStamps_Tn = tStamp(tn_idx);
            if(TU)
                TimeStamps_TU = tStamp(~(tp_idx | tn_idx)); 
            end
        % Iterate to find 50 random response variables using M1 indirects
        for iter=(1:50)
            clear tp_unit tn_unit           
            display(iter);

            if ~TU

                TimeStamps_Tp_Iter = TimeStamps_Tp(randperm(length(TimeStamps_Tp), tp_count));
                TimeStamps_Tp_Iter = reshape(TimeStamps_Tp_Iter,1,[]);
                tp_idx = 1:tp_count;
                tp_unit = fn_getbinnedData(TimeStamps_Tp_Iter,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,1,tp_idx);   

                if tn_count              
                    TimeStamps_Tn_Iter = TimeStamps_Tn(randperm(length(TimeStamps_Tn), tn_count));
                    TimeStamps_Tn_Iter = reshape(TimeStamps_Tn_Iter,1,[]);
                    tn_idx = 1:tn_count;
                    tn_unit = fn_getbinnedData(TimeStamps_Tn_Iter,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,1,tn_idx);
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
                res{b,n,iter} = fn_getbinnedData(TimeStamps_TU_Iter,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,1,tu_idx);  
            end
 
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
        id_idx1 = (3*std_id + mean_id < max_id);
%         id_idx2 = ( mean_d - 2.5*std_id > min_id);
%         idx = id_idx1 | id_idx2;
        idx = id_idx1;
        
        
        indirect_unit = fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,trial_start,trial_stop,bin_window,0,idx);
                  
        indirect{b,n} = indirect_unit;
        
        
        clear all_trial* TimeStamps2 rewards_onse*
%         end
   end   

end

%% Train the GLM model and test


  all_R2 = zeros(length(bin),50,20);
  r2 = zeros(length(bin),50,20);

rootpath = 'Z:\Data_for_SpikeGLM\Results';

  fileName = '\GLM_CbInd_M1TU_TaskStart_allData.mat';
  m = matfile([rootpath,fileName],'Writable',true);
  m.r2 = r2;
  m.all_R2 = all_R2;
%   m.all_model = all_model;
  
%iterate for all binned samples  
for b=2:4%1:4%length(bin)    
    
% iterate for 10 time for crossvalidation    
for iter=1:50

    disp(iter);
    res_cell  =  res(b,first:last,iter);
    indirect_cell  =  indirect(b,first:last);
  
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
          for t=floor(dur/2):size(res_data,2)-(floor(dur/2)+1)
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
%           m.all_model(b,iter,i) = mdl;
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

 save([rootpath,fileName],'all_model','-append');

%% Compare the R2 for Cb ind - M1 ind and Cb ind - M1 dir


load([rootpath,'\GLM_model_result_M1d_CbId_AroundReward_new.mat']);
    avg_1 = mean(all_R2,2);
    clear all_R2;
load([rootpath,'\GLM_CbInd_M1Ind_results_AroundReward.mat']);
    avg_2 = mean((all_R2),2);
for j=1:20%size(m.all_R2,3)
    subplot(10,2,j);
    bar([avg_1(4,:,j),avg_2(4,:,j)]);
end
%% Compare the R2 for Cb ind - M1 ind and Cb ind - M1 dir

load([rootpath,'\GLM_model_result_M1d_CbId_AroundReward.mat']);
    avg_1 = mean(all_R2,2);
    clear all_R2;
load([rootpath,'\GLM_CbInd_M1Ind_results.mat']);
    avg_2 = mean((all_R2),2);
for j=1:20%size(m.all_R2,3)
    subplot(10,2,j);
    bar([avg_1(4,:,j),avg_2(4,:,j)]);
end

%% Compare the R2 for Cb ind - M1 ind and Cb ind - M1 TU
rootpath = 'Z:\Data_for_SpikeGLM\Results';
load([rootpath,'\GLM_model_result_M1d_CbId_TaskStart_allData.mat'],'all_R2');
    avg_1 = mean(all_R2,2);
    clear all_R2;
load([rootpath,'\GLM_CbInd_M1Ind_TaskStart_allData.mat'],'all_R2');
    avg_2 = mean(all_R2,2);
    clear all_R2;
load([rootpath,'\GLM_CbInd_M1TU_TaskStart_allData.mat'],'all_R2');
    avg_3 = mean((all_R2),2);
for j=1:20%size(m.all_R2,3)
    subplot(10,2,j);
%     bar([avg_1(4,:,j),avg_2(4,:,j),avg_3(4,:,j)]);
    bar([avg_1(3,:,j),avg_2(3,:,j),avg_3(3,:,j)]);
end

figure;
% plot([1:length(avg_1)],avg_1(3,:),'-o'); hold on;
% plot([1:length(avg_2)],avg_2(3,:),'-x'); hold on;
% plot([1:length(avg_3)],avg_3(3,:),'-s'); hold on;
% 
plot([1:length(avg_1)],avg_1(4,:),'-o'); hold on;
plot([1:length(avg_2)],avg_2(4,:),'-x'); hold on;
plot([1:length(avg_3)],avg_3(4,:),'-s'); hold on;
legend(["Cb TRi - M1 TRd","Cb TRi - M1 TRi","Cb TRi - M1 TU"]);
xlabel("Session");
ylabel("R2");

%% Create mat file with R2 data

TRd_R2 = avg_1(4,:);
TRi_R2 = avg_2(4,:);
TU_R2 = avg_3(4,:);
save([rootpath,'\R2_data.mat'], 'TRd_R2', 'TRi_R2', 'TU_R2');
display("saved");

%%
load([rootpath,'\R2_data.mat'], 'TRd_R2', 'TRi_R2', 'TU_R2');
figure;
cons = [TRd_R2;TRi_R2;TU_R2]
boxplot(cons');
%% Mixed effect model
Y = [avg_1(4,:), avg_2(4,:), avg_3(4,:)];
mdl = [ones(1,20),ones(1,20)*2,ones(1,20)*3];
mdl = double(Unit);
Session =[(1:20),(1:20),(1:20)];
Session = double(Session);
lme1 = fitlmematrix(mdl',Y',ones(size(Y')),Session')
lme2 = fitlmematrix(mdl',Y',Session',[])
%%
figure;
cons = [(avg_1(4,:));(avg_2(4,:));(avg_3(4,:))];

bar(nanmean(cons,2)); hold on;
errorbar([1,2,3],nanmean(cons,2),nanstd(cons,[],2)/sqrt((size(cons,2)-1))); hold on;
plot(cons,'k')
xticklabels({"Cb TRi - M1 TRd","Cb TRi - M1 TRi","Cb TRi - M1 TU"});
%%
figure;
cons = [(avg_1(4,:));(avg_2(4,:));(avg_3(4,:))];

boxplot((cons')); hold on;
% errorbar([1,2,3],nanmean(cons,2),nanstd(cons,[],2)/sqrt((size(cons,2)-1))); hold on;
xticklabels({"Cb TRi - M1 TRd","Cb TRi - M1 TRi","Cb TRi - M1 TU"});
%%
figure;
% cons = [(avg_1(3,:));(avg_2(3,:));(avg_3(3,:))];
cons = [(avg_1(4,:));(avg_2(4,:));(avg_3(4,:))];
len = length(cons);
x = [ones(1,len),2*ones(1,len),3*ones(1,len)];
bar(nanmean(cons,2)); hold on;
errorbar([1,2,3],nanmean(cons,2),nanstd(cons,[],2)/sqrt((size(cons,2)-1))); hold on;
cons = cons';
scatter(x,cons(:));hold on;
ylabel("R2");
xticklabels({"Cb TRi - M1 TRd","Cb TRi - M1 TRi","Cb TRi - M1 TU"});

%%
figure;
violinplot(cons');
ylabel("R2");
xticklabels({"Cb TRi - M1 TRd","Cb TRi - M1 TRi","Cb TRi - M1 TU"});
%%
%%

figure;
plot([1:length(avg_1)],avg_1(3,:),'-o'); hold on;
plot([1:length(avg_2)],avg_2(3,:),'-x'); hold on;
legend(["Cb TRi - M1 TRd","Cb TRi - M1 TRi"]);
xlabel("Session");
ylabel("R2");
%%
mean(avg_1(4,:))

mean(avg_2(4,:))
%%
avg_R2= mean(abs(all_R2),2);
figure;
bar((avg_R2(:,1,:)));
figure;
bar(mean(abs(r2),2));

%% Plot all R2
    figure;
    for j=1:size(all_R2,3)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
%         tmp = cell2mat(r2);
        tmp = avg;
        if(sum(tmp))
            plot(tmp(1:4),'.-','MarkerSize',15,'LineWidth',2); hold on;
%             count = count +1;
        end
    end
%%
figure;
model = all_model{5,6};
bar(model.Coefficients.Estimate(2:end));
%% Get indirect weights and probability for all sessions - M1 BMI 

% for [j,k]=[1:10,1:10]
rootpath = 'Z:\Data_for_SpikeGLM\Results\';
% fileName = 'GLM_CbInd_M1Ind_TaskStart_allData.mat';
av = mean(all_R2(4,:,:),2);
fileName = 'GLM_model_result_M1d_CbId_TaskStart_allData.mat';
%     
    clear all_model all_R2 r2;
    load([rootpath, fileName]);
%     for j=18:19%1:size(all_model,2)
    for j=10%19%1:size(all_model,3)
    all_bins = [];
%   for k=2%1:size(all_model,1)
    for k=2%1:size(all_model,2)
        if(av(j)>0)
%         if(all_R2(4,k,j)>0)
        figure;
%         model = all_model{k,j};
        model = all_model{2,k,j};
        weights = model.Coefficients.Estimate(2:end);

        w = reshape(weights,length(weights)/11,11);
%         weights = reshape(abs(weights),length(weights)/17,17);
        weights = reshape(abs(weights),length(weights)/11,11);
        weights = weights(any(weights,2),:);

  
        weight = normalize(weights,2,'norm',inf);

        [m,idx] = nanmax(weight,[],2);
        [ii,idx] = sort(idx); 

        weight = weight(idx,:);
        % weight = sortrows(weight,max(weight'));
        subplot(2,1,1);
        image(weight, 'CDataMapping', 'scaled'); hold on;
        colorbar; hold on;
        caxis([0,1]);
        xticks((1:11));

        xticklabels(["-50","-40","-30","-20",'-10','0','10','20','30','40','50']);
        xlabel("t(ms)");
        ylabel("Cb Units");
%         set(gcf,'WindowState','maximized')

%         saveas(gcf,[rootpath,'\Fig\','IndirectWeights_M1TRd_Session_',num2str(j),'_',num2str(k),'.fig']);
%         close;
        all_bins{k}=ii;
        end
    end
    
    [GC,GR] = groupcounts(cell2mat(all_bins(:))); 
%     histogram(cell2mat(all_bins(:)));
    if(~isempty(GC))
        
%     figure;
    subplot(2,1,2);
    bar(GR,GC/1.0,1);
    xlim([0.5,11.5]);
    xticks((1:11));
    xticklabels(["-50","-40","-30","-20",'-10','0','10','20','30','40','50']);
    xlabel("t(ms)");
    ylabel("Count");
%     xticklabels(['-50','-40','-30','-30','-20','-10','0','10','20','30','40','50']);

    end
    end

%% Create mat file with weights data
TRd_weights = weights;
TRd_counts = GC;

save([rootpath,'\TRd_GLM_data.mat'], 'TRd_weights', 'TRd_counts');

display("saved");
%% Create mat file with weights example
w1 = w(25,:);
w2 = w(50,:);
save([rootpath,'\TRd_example_weights.mat'],'w1','w2');
%% Plot weights

figure;
subplot(2,1,1);
bar(w(25,:));
xticks((1:11));
xticklabels(["-50","-40","-30","-20",'-10','0','10','20','30','40','50']);
xlabel("t(ms)");
ylabel("Weight");        
% yline(0);
subplot(2,1,2);
bar(w(50,:));
xticks((1:11));
xticklabels(["-50","-40","-30","-20",'-10','0','10','20','30','40','50']);
xlabel("t(ms)");
ylabel("Weight");
% yline(0);

%% Get highest weight time lag histogram - M1 BMI

bin = [5,10,25,50,100];
% for [j,k]=[1:10,1:10]

% for [j,k]=[1:10,1:10]
rootpath = 'Z:\Data_for_SpikeGLM\Results\';
% fileName = 'GLM_CbInd_M1Ind_TaskStart_allData.mat';
fileName = 'GLM_model_result_M1d_CbId_TaskStart_allData.mat';
all_bins = [];

%     clear all_model all_R2;
%     load([rootpath, fileName]);
    for j=1:size(all_model,3)
%     for j=1:size(all_model,2)        
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(4)~=0)
            
        for k=1:10
            for b=3%1:4%2%2:6
%                 figure;
                model = all_model{b,k,j};
%                 model = all_model{b,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                all_bins{j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
        end
        end
    end

%%
    figure;
    
    [GC,GR] = groupcounts(cell2mat(all_bins(:))); 
%     histogram(cell2mat(all_bins(:)));
    bar(GC/10.0);
%     saveas(gcf,[rootpath,'\','Histogram_stroke.tif']);
%     close;
%%    
 figure
 pd = fitdist(cell2mat(all_bins(:)),'Kernel','Kernel','epanechnikov');
% pd = histfit(cell2mat(all_bins(:)));
xlim([0,18]);
x = 0:.1:18;
ySix = pdf(pd,x);
plot(x,ySix,'k-','LineWidth',2)
    
%%
    figure;
    

    bin_2 = all_bins(:,:,3);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    yline(mean(counts)); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;
