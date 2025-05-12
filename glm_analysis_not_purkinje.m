% GLM model
% Predictors were binned firing rates of either all M2 units or all M1 indirect units


%% Generate binned data

clc;clear;

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
    [0,0,0, 9, 6, 11, 4, 45,19, 11, 6, 1, 20,43, 17,23],...
    [0,0,0,44, 27, 57, 20, 4, 1,  33, 1 ]};

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
    [0,0,0,85,80,93,100,82,140,55,100,38,60,120,60,82],...
    [0,0,0,62, 83, 75, 48, 61,141,113,35]};



sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160', 'I161'};


stroke= [4:5,9:13];

intact = [1:13];
intact = setdiff(intact, stroke);
Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';


robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6],[8,14,15,16,18,19],[7,9,11,14],[8,9,10]};




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
for b=2:4%length(bin)
for s=12%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    load([rootpath,'\purkinje_units_161.mat'])
    load([rootpath,'\purkinje_ch_161.mat'])
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1
    
        if isfield(channels,bmiBlocks(n).name)
    % Define save path
        savepath = [rootpath, sub{s},'\Data\',bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');
        events = rewards_onset;
        
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

        bin_window = bin(b)/1000; % in seconds
%         tp = [tp;fn_getbinnedData(TimeStamps_tp,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
%         tn = [tn;fn_getbinnedData(TimeStamps_tn,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
        

        
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
                 % Not Purkinje putative
                 idx = find(ch==k+1);
                 if ~isempty(idx) && unit(idx(1)) ~= 0                    
                 else
                    TimeStamps(j,k) = TimeStamps2(j,k);
                 end   
             end
          end
        end
        
%       indirect = [indirect;fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart(n),tstop(n),bin_window)];
        indirect_unit = fn_getbinnedData(TimeStamps,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
                  
        indirect{b,n} = indirect_unit;
        
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
   end   

end
end

% Train the GLM model and test

%   lag = [-200,-180,-160,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200];
%   
  all_R2 = zeros(length(bin),10,10);
  
  r2 = zeros(length(bin),10,10);
  prob = []; 
  all_model = [];
  
  
  
%iterate for all binned samples  
for b=2:4%length(bin)    
    
% iterate for 10 time for crossvalidation    

for iter=1:10
    display(iter);
%       
%     tp_cell  =  tp(b,:);
%     tp_cell = tp_cell(~cellfun('isempty',tp_cell));
%     tn_cell  =  tn(b,:);
%     tn_cell = tn_cell(~cellfun('isempty',tn_cell));

    res_cell  =  res(b,:);
    res_cell = res_cell(~cellfun('isempty',res_cell));
%     indirect_data  =  reshape(indirect{b},1,size(indirect{b},1)*size(indirect{b},2));
    indirect_cell  =  indirect(b,:);
    indirect_cell = indirect_cell(~cellfun('isempty',indirect_cell));
%     indirect_data = cell2mat(indirect_data(:));
    indirect_data = [];
%     tp_data = [];
%     tn_data = [];
    res_data = [];
    
    
    
    for i=1:length(indirect_cell)
        indirect_data = [];

        res_data = [];

        indirect_data1 = cell2mat(indirect_cell{i}(:));
        res_data = cell2mat(res_cell(i));
        num_trial = size(res_data, 1);    
        
%         load([rootpath, sub{s},'\Data\', bmiBlocks(first(s)+i-1).name, '\SigInd.mat'], 'idx');
%         idx_ = repmat(idx,num_trial,1);
%         indirect_data =  indirect_data1(idx_,:);
        
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
%             Xnew = pred((8*end/10)+1:end,:);
%             z = var((8*end/10)+1:end,:);
%             Xnew = pred(:,floor(end/10)+1:end); 
            ypred = predict(mdl,Xnew);
            
            z = var((9*end/10)+1:end,:);
            r = z-ypred;
            normr = norm(r);
            SSE = normr.^2; 
            SST = norm(z-mean(z))^2;
%             sumr = sum(r);
%             SSE = sum((z-ypred).^2); 
%             SST = sum((z-mean(z)).^2);
%             SST = (length(z)-1)*var(z);
            R2 = 1 - SSE/SST;
%             for i=1:9
%               Xnew = pred(:,i*end/10+1:((i+1)*end)/10);
%               ypred = predict(mdl,Xnew);
%     %           y{i} = ypred;
% 
%                 z = var;
%                 z_est = ypred;
%                 r = z-z_est;
%                 normr = norm(r);
%                 SSE = normr.^2;
%                 SST = norm(z-mean(z))^2;
%                 R2{i+1} = 1 - SSE/SST;
%             end
            all_R2(b,iter,i) = R2;
        end
      end
end

end

% Save the model/result and plot

save([rootpath, sub{s}, '\GLM_model_result_not_purkinje.mat'], 'all_model','r2','all_R2', '-v7.3');
display("saved");

% PLot R2 and save
figure;
for j=1:size(all_R2,3)
    avg = mean((all_R2),2);
    subplot(5,2,j);
    bar(avg(:,:,j));
end


%%
avg_R2= mean(abs(all_R2),2);
figure;
bar((avg_R2));
figure;
bar(mean(abs(r2),2));
%%
figure;
model = all_model{5,6};
bar(model.Coefficients.Estimate(2:end));


%% Get indirect weights and probability for all sessions

% for [j,k]=[1:10,1:10]
rootpath = 'Z:\Rohit\BMI_Data\';

sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154'};
for s=11 
    
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_result_purkinje.mat']);
    for j=1:size(all_model,3)
    all_bins = [];
    for k=1:10
        figure;
        model = all_model{2,k,j};
        weights = model.Coefficients.Estimate(2:end);
%         weights = reshape(abs(weights),length(weights)/17,17);
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

        saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_new_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
        close;
        all_bins{k}=ii;
    end
    figure;
    [GC,GR] = groupcounts(cell2mat(all_bins(:))); 
%     histogram(cell2mat(all_bins(:)));
    bar(GC/10.0);
    saveas(gcf,[rootpath,'\',sub{s},'\','Histogram_new_Sesssion_',num2str(j),'.tiff']);
    close;
    end


    for j=1:size(all_R2,3)
        avg = mean((all_R2),2);
        subplot(5,2,j);
        bar(avg(:,:,j));
    end
    saveas(gcf,[rootpath,'\',sub{s},'\','Test_R2_purkinje.tiff']);
    close;
    figure;
    for j=1:size(all_R2,3)
        avg = mean((r2),2);
        subplot(5,2,j);
        bar(avg(:,:,j));
    end
    saveas(gcf,[rootpath,'\',sub{s},'\','Training_R2_purkinje.tiff']);
    close;
end


%% Get highest weight time lag histogram

bin = [5,10,25,50,100];
% for [j,k]=[1:10,1:10]
rootpath = 'Z:\Rohit\BMI_Data\';
% sub = {'I096','I107','I110','I111','I112'};
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160','I161'};
all_bins = [];
stroke = [4:5,9:13];
for s=intact
    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result_not_purkinje.mat']);
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
                all_bins{s,j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
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
    

    bin_2 = all_bins(:,:,:,2);
    counts = histcounts(cell2mat(bin_2(:)),'Normalization','probability');
    line([0,12],[mean(counts),mean(counts)]); hold on;
    ylim([0,0.25]);
    bar(counts); hold on;
%     for i=1:4
% 
%         temp = all_bins(:,:,:,i);
%         tmp = max(cell2mat(temp(:)));
%         pd = fitdist(cell2mat(temp(:)).*(11/(tmp+1)),'Kernel');
% 
%         % pd = histfit(cell2mat(all_bins(:)));
%         xlim([0,12]);
%         x = 0:.05:12;
%         ySix = pdf(pd,x);
%         [f,xi] = ksdensity(ySix,'Bandwidth',0.1);
% %         ySix = ySix.*(12/(tmp+1));
% %         plot(x,ySix,'k-','LineWidth',2); hold on;
%         plot(xi,f,'k-','LineWidth',2); hold on;
% %         pbaspect([(12/(tmp+1)) 1])
%     end;

%     saveas(gcf,[rootpath,'\','Histogram_prob_healthy.tiff']);
%     close;
% 
%     for j=1:size(all_R2,3)
%         avg = mean((all_R2),2);
%         subplot(5,2,j);
%         bar(avg(:,:,j));
%     end
%     saveas(gcf,[rootpath,'\',sub{s},'\','Test_R2.tiff']);
%     close;
%     figure;
%     for j=1:size(all_R2,3)
%         avg = mean((r2),2);
%         subplot(5,2,j);
%         bar(avg(:,:,j));
%     end
%     saveas(gcf,[rootpath,'\',sub{s},'\','Training_R2.tiff']);
%     close;

%% R2 healthy vs stroke
stroke= [4:5,9:13];

intact = [1:13];
intact = setdiff(intact, stroke);
rootpath = 'Z:\Rohit\BMI_Data\';
robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6],[8,14,15,16,18,19],};

sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128', 'I154', 'I160','I161'};
% 4 sessions are not significant glm sessions
% sig_session = {[4,5],[3,4,6],[7,8,9,10,13,14],[3,4,5,8,9,10],[5,7,9],...
%     [4,6,8],[5:9],[8],[8],[3,6],[8,14,15,16,18]};
sig_session = {[1,1,0,0,0,0,0,0,0,0],[1,1,1,0,0,0,0,0,0,0],[0,1,1,1,1,1,1,0,0,0],[1,1,1,1,1,1,0,0,0,0],[1,0,1,1,0,0,0,0,0,0],...
    [1,1,1,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[1,0,1,0,0,0,0,0,0,0],[1,1,1,1,1,0,0,0,0,0]};

val_healthy = [];
val_stroke = [];
for s=intact
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_result_not_purkinje.mat'],'all_R2');
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
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_result_not_purkinje.mat'],'all_R2');
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
bar([mean(val1),mean(val2)]);
errorbar([1,2],[mean(val1),mean(val2)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1))]);
scatter([repelem(1,length(val1)),repelem(2,length(val2))],[val1;val2]);
[h,p] = ttest2(val1,val2);

%%
group = [repelem(1,length(val1)),repelem(2,length(val2))];
[p,t,stats] = anovan([val1;val2]',{group});



% [c,m,h,gnames] = multcompare(stats); 
% saveas(gcf,[rootpath,'\','R2-healthyvsStroke_new.tiff']);
%     close;

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
formula = 'R2 ~ 1 + IS + (1 | RatID)';
lme_R2 = fitlme(tbl,formula)

%% PLot all R2 for all session

rootpath = 'Z:\Rohit\BMI_Data\';
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154'};
robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6],[8,14,15,16,18,19]};
val_healthy = [];
val_stroke = [];
r2  = {};
figure;
xlim([0,4]); hold on;
xticks([0,1,2,3,4]);
xticklabels({'','10 ms','25 ms', '50 ms',''});
ylabel('R2');
count = 0;
for s=stroke

    clear all_model all_R2;
    load([rootpath, sub{s}, '\GLM_model_result_new.mat'],'all_R2');
%     for j=1:size(all_R2,3)
        
    for j=1:length(robust_session{s})
        avg = mean((all_R2),2);
        r2{s,j} = avg(:,:,j);
        tmp = cell2mat(r2(s,j));
        if(sum(tmp))
            plot(tmp(2:4),'.-','MarkerSize',15,'LineWidth',2); hold on;
            count = count +1;
        end
    end
%     if(~isempty(r2{s,j}))
%         scatter(r2{s,j}); hold on;
%     end
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

%% PLot R2 and save
figure;
for j=1:size(all_R2,3)
    avg = mean((all_R2),2);
    subplot(5,2,j);
    bar(avg(:,:,j));
end
%%
% saveas(gcf,[rootpath,'\',sub{s},'\','Test_R2.tiff']);
% close;
figure;
for j=1:size(all_R2,3)
    avg = mean((r2),2);
    subplot(5,2,j);
    bar(avg(:,:,j));
end
figure;
for j=1:size(all_R2,3)
    avg = mean((adj_r2),2);
    subplot(5,2,j);
    bar(avg(:,:,j));
end
% saveas(gcf,[rootpath,'\',sub{s},'\','Training_R2.tiff']);
% close;
% set colormalt to "jet"
% cmap = jet(10);
% colormap(cmap)

% bar(model.Coefficients.Estimate(2:end));

%%
figure;
histogram(ii)
%%
for k=1:length(all_R2)
    figure;
    bar(abs(cell2mat(all_R2{k}))); 
end


%%

avg_R2 = [];
for i=1:length(all_R2)
    avg_R2{i} = mean(abs(cell2mat(all_R2{i}))); 
%         avg_R2{i} = mean((cell2mat(all_R2{i}))); 

end
bar(cell2mat(avg_R2));

%%

%
% save([rootpath,'\GLM_model_data.mat'], 'indirect','tp','tn','all_R2','r2');

  %%

z = res';
z_est = ypred;
r = z-z_est;
normr = norm(r);
SSE = normr.^2;
SST = norm(z-mean(z))^2;
R2 = 1 - SSE/SST;

%%


p = coefTest(mdl)
tbl = devianceTest(mdl)
%plotSlice(mdl);

%%
% confint = coefCI(mdl)
% plotDiagnostics(mdl)

plotResiduals(mdl)
figure;
plotResiduals(mdl,'fitted')
figure;
plotResiduals(mdl,'probability')


%% Potent space (Cb direct)
clear; clc; close all;
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\I111\Data\';
bmiBlocks = dir(rootpath);
before_zero = 1.5;
after_zero =0.5;


%%


%     indirect_data  =  reshape(indirect{b},1,size(indirect{b},1)*size(indirect{b},2));
    indirect_data  =  indirect(b,:);
indirect_data = indirect_data(~cellfun('isempty',indirect_data));
% indirect_data = cell2mat(indirect_data(:));
A = [];
for i=1:length(indirect_data)
    A = [A;cell2mat(indirect_data{i}(:))];
end





%   fitglm()