% GLM model - for M1 predicting purkinje - GLM-1p
% Predictors were binned firing rates of M1 indirects
% Rohit 2023


%% Generate binned data

clc;clear;

%Load all session info
bmi_session_info;

% bmiBlocks = dir(rootpath); 
before_zero = 2.2;
after_zero =0.2;

tp = [];
tn = [];

res = [];
indirect = [];
    
bin = [5,10,25,50,100]; % binnning at different ms
TU = 0;

for z=11:15 

bin = [5,10,25,50,100]; % binnning at different ms
for b=2:4%1:4%length(bin)
for s=z%M1
    s
    rootpath = 'Z:\Rohit\BMI_Data\';
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    load([rootpath,'\purkinje_units_all.mat'])
    load([rootpath,'\purkinje_ch_all.mat'])
    for n=intersection{s}
    
        
    % Define save path
        savepath = [rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];
        load([rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Timestamps_Direct.mat']);
        load([rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');

        events = rewards_onset;
        load([rootpath,sub{s},'\Data\' bmiBlocks(n).name,'\Performance_stats_early_late.mat']);

        bin_window = bin(b)/1000; % in seconds

        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);

        
        % Get only good units
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
        for iter=(1:50)
            clear tp_unit tn_unit           
            display(iter);


            if ~isempty(TimeStamps_Tp)

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
   end   
end
end

% Train the GLM model and test


  all_R2 = zeros(length(bin),50,20);
  r2 = zeros(length(bin),50,20);

  fileName = '\GLM_1P.mat';
  m = matfile([rootpath, sub{s},fileName],'Writable',true);
  m.r2 = r2;
  m.all_R2 = all_R2;
%   m.all_model = all_model;
  
%iterate for all binned samples  
for b=2:4%1:4%length(bin)    
    
% iterate for 10 time for crossvalidation    

for iter=1:50

    disp(iter);

    res_cell  =  res(b,intersection{s},iter);
    indirect_cell  =  indirect(b,intersection{s});
  
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
%           mdl = fitglm(pred(1:(8*end/10),:), var(1:(8*end/10)));  %% response variable is the data at reward
%           m.all_model(b,iter,i) = mdl;
          if(iter==1)
              all_model{b,i} = mdl;
          end
          %R-squared value for model
          m.r2(b,iter,i) = mdl.Rsquared.Ordinary;
%           adj_r2(b,iter,i) =  mdl.Rsquared.Adjusted;
%           prob{b,iter,i} = mdl.Fitted.Response;
          
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

%% Get highest weight time lag histogram

%Load all session info
bmi_session_info;
str_bins = [];
int_bins = [];
for s=M1_stroke
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat']);
    for j=1:size(all_model,2)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
%         for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                str_bins{s,j,b}=ii;
%                 tot_session = tot_session + 1;
            end    
%         end
        end
    end
end 



for s=[1,2,6,7]%M1_intact
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat']);
    for j=1:size(all_model,2)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
%         for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                int_bins{s,j,b}=ii;
%                 tot_session = tot_session + 1;
            end    
%         end
        end
    end
end 

%% plot histogram with kernel density estimates - intact
    figure;
    
    bin_2 = int_bins(:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    yline(mean(counts)); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;
     for i=2:4

        temp = int_bins(:,:,i);
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
    
    bin_2 = str_bins(:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    yline(mean(counts)); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;
     for i=2:4

        temp = str_bins(:,:,i);
        tmp = max(cell2mat(temp(:)));
        pd = fitdist(cell2mat(temp(:)).*(12/(tmp+1)),'Kernel');

        xlim([0,12]);
        x = 0.5:.05:12;
        ySix = pdf(pd,x);
        plot(x,ySix); hold on;
    end
    legend('uniform', 'histogram', '10 ms', '25 ms', '50 ms');

 %% ks test
 str_bin = str_bins(:,:,2);
 str_counts = histcounts(cell2mat(str_bin(:)),'Normalization','probability');
 int_bin = int_bins(:,:,2);
 int_counts = histcounts(cell2mat(int_bin(:)),'Normalization','probability');
 
%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(str_bin(:)), cell2mat(int_bin(:)), 1000)
 %%
 bar(int_counts); hold on;
 yline(mean(int_counts)); hold on;
 figure;
 bar(str_counts); hold on;
 yline(mean(int_counts)); hold on;
 
 %% check for uniformity

[h,p,ksstat] = kstest2(int_counts,repelem(mean(int_counts),11))
[h,p,ksstat] = kstest2(str_counts,repelem(mean(str_counts),11))

% [h,p] = chi2gof(int_counts',[1:11])
% [h,p] = chi2gof(str_counts',[1:11])
%% bootstrap

[ci,bootstat] = bootci(100000,@immse,int_counts,repelem(mean(int_counts),11));
ci
hist(bootstat); hold on;
vline(ci(2)); hold on;
vline(immse(int_counts,repelem(mean(int_counts),11)),'k');
[h,p]= ttest(bootstat)
figure;
[ci,bootstat] = bootci(100000,@immse,str_counts,repelem(mean(str_counts),11));
ci
hist(bootstat); hold on;
vline(ci(2)); hold on;
vline(immse(str_counts,repelem(mean(str_counts),11)),'k');
[h,p]= ttest(bootstat)
%%
    figure;
    

    bin_2 = all_bins(:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    yline(mean(counts)); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;
    

%% R2 healthy vs stroke
bmi_session_info;

val_healthy = [];
val_stroke = [];
for s=[1,2,6,7]
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat'],'all_R2');
    for j=1:size(all_R2,3)
%     if (sig_session{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        
        if(~isnan(r2(4)) && r2(4))
            val_healthy{s,j} = r2(4);    
        end
%     end
    end
end 

for s=M1_stroke
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat'],'all_R2');
    for j=1:size(all_R2,3)
%     if (sig_session{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(~isnan(r2(4)) && r2(4) )
            val_stroke{s,j} = r2(4);    
        end
%     end
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

