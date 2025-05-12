% GLM model - for purkinje CB direct - GLM-pd
% Predictors were binned firing rates of purkinje
% Rohit 2023

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
for b=2:4%length(bin)
for s=15%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    load([rootpath,'\purkinje_units_all.mat'])
    load([rootpath,'\purkinje_ch_all.mat'])
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1
    
        if isfield(channels,bmiBlocks(n).name)
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
        
        indirect_unit = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);         
        indirect{b,n} = indirect_unit;
 
        end
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
for b=2:4%length(bin)    
    
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
       indirect_data =  indirect_data1;
        
        if ~isempty(indirect_data)
        num_indirect =floor(size(indirect_data, 1)/size(res_data, 1));

          start = 0;
          
%           dur = 400/bin(b);
          dur = 100/bin(b);
%           stop = dur -1; 
          stop = dur; 
          pred = [];
          var = [];
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

% Save the model/result and plot

save([rootpath, sub{s}, '\GLM_model_result_purkinje.mat'], 'all_model','r2','all_R2', '-v7.3');
display("saved");

% PLot R2 and save
figure;
for j=1:size(all_R2,3)
    avg = mean((all_R2),2);
    subplot(5,2,j);
    bar(avg(:,:,j));
end


%% Get highest weight time lag histogram

bin = [5,10,25,50,100];

bmi_session_info;

for s=stroke
    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result_purkinje.mat']);
    for j=1:size(all_model,3)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
            
        for k=1:10
            for b=2:4%2%2:6
%                 figure;
                model = all_model{b,k,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

                str_bins{s,j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
        end
        end
    end
end 


for s=intact
    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result_purkinje.mat']);
    for j=1:size(all_model,3)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
            
        for k=1:10
            for b=2:4%2%2:6
%                 figure;
                model = all_model{b,k,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
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
    ylim([0,0.5]);
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
    ylim([0,0.5]);
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
 %% ks test
 str_bin = str_bins(:,:,:,2);
 str_counts = histcounts(cell2mat(str_bin(:)),'Normalization','probability');
 int_bin = int_bins(:,:,:,2);
 int_counts = histcounts(cell2mat(int_bin(:)),'Normalization','probability');
 
%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(str_bin(:)), cell2mat(int_bin(:)), 10000)
 
 bar(int_counts); hold on;
 bar(str_counts); hold on;
 
 %% check for uniformity

[h,p,ksstat] = kstest2(int_counts,repelem(mean(int_counts),11))
[h,p,ksstat] = kstest2(str_counts,repelem(mean(str_counts),11))


% bootstat = bootstrp(10000,@immmse, int_counts,repelem(mean(int_counts),11));
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

    bin_2 = all_bins(:,:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    yline(mean(counts)); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;

%% R2 healthy vs stroke
bmi_session_info;

val_healthy = [];
val_stroke = [];
for s=intact
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_result_purkinje.mat'],'all_R2');
    for j=1:size(all_R2,3)
%     if (sig_session{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        
        if(r2(4))
            val_healthy{s,j} = r2(4);    
        end
%     end
    end
end 

for s=stroke
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_result_purkinje.mat'],'all_R2');
    for j=1:size(all_R2,3)
%     if (sig_session{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(4))
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
errorbar([1,2],[mean(val1),mean(val2)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1))],'k','linewidth',4);
scatter([repelem(1,length(val1)),repelem(2,length(val2))],[val1;val2],'k');
[h,p] = ttest2(val1,val2)
set(gca,'TickDir','out');
xlim([0.5,2.5]);
ylim([0,0.7]);

%%
group = [repelem(1,length(val1)),repelem(2,length(val2))];
[p,t,stats] = anovan([val1;val2]',{group});

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
