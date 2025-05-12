%% Get spontaneuos activity spectrum data for stroke animals
%  Author: Rohit Rangwani - 2025
%  ---------------------------------------------------------------------
%% Read collated LFP file for Cb and analyze

clear;clc;close all; tic;

bmi_session_info;

for s=intact
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Collated_LFP_new1.mat'],'file')

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Collated_LFP_new1.mat'],'badchans2');
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trials');

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_Cb.mat']);


       % Cb
        data = lfp_Cb;
        data(badchans2,:)  = [];

        event = all_trials(1);
        % get 2 min spontaneous activity before session start
        if(event-60*Fs > 0)
            data = data(:,event-60*Fs:event);
        end

        % Clear variables 
      clear trial_data_n ersp_data trial_data1_erp trial_data1_cmr
      
      % Z-scoring
      for j=1:size(data,1)
        tmp = reshape(squeeze(data(j,:)),1,[]);
        trial_data_n(j,:,:) = (data(j,:,:)-mean(tmp))./std(tmp);
      end
      
      % Median subtraction
      med = median(trial_data_n,1);
      trial_data1_cmr = bsxfun(@minus,trial_data_n,med);
      
      % Get Event=Related Spectral Perturbation (ERSP) using eeglab for M1 without ERP subtraction
      for j=1:size(trial_data1_cmr,1)
        disp(['Ch_',num2str(j)]);
        data = squeeze(trial_data1_cmr(j,:,:));
        [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,size(data,1),[0 size(data,1)-1],Fs,[0.01 0.1],...
          'baseline',NaN,'freqs', [0 10],'verbose','off','trialbase','on','plotersp','off','plotitc','off');
        clf;
      end
      
      % Get normalized ESRP data generated before subtracting ERP
      data     = ersp_data.*conj(ersp_data);
      bl       = mean(data,3);
      bl_std   = std(data,[],3);
      datanorm = bsxfun(@minus,data,bl);
      ersp_datanorm1_cmr = bsxfun(@rdivide,datanorm,bl_std);
      ersp_data1_cmr = data;
      
      % ERP subtraction (Channel by Channel)
      for ch=1:size(trial_data1_cmr,1)
        single_channel_data = squeeze(trial_data1_cmr(ch,:,:));
        mean_erp = mean(single_channel_data,1);
        trial_data1_erp(ch,:,:) = bsxfun(@minus,single_channel_data,mean_erp);
      end
      
      % Get Event=Related Spectral Perturbation (ERSP) using eeglab for M1 with ERP subtraction
      for j=1:size(trial_data1_erp,1)
        disp(['Ch_',num2str(j)]);
        data = squeeze(trial_data1_erp(j,:,:));
        [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,length(data),[1 length(data)],Fs,[0.01 0.1],...
          'baseline',NaN,'freqs', [0 10],'verbose','off','trialbase','on','plotersp','off','plotitc','off');
        clf;
      end
      %,'plotersp','off','plotitc','off'
      % Get normalized ESRP data generated after subtracting ERP
      data = ersp_data.*conj(ersp_data);
      bl       = mean(data,3);
      bl_std   = std(data,[],3);
      datanorm = bsxfun(@minus,data,bl);
      ersp_datanorm1_erp = bsxfun(@rdivide,datanorm,bl_std);
      ersp_data1_erp = data;

      delta_Cb = squeeze(ersp_data1_erp(:,freqs>=0.4&freqs<=4,:,:));
      delta_Cb = squeeze(mean(mean(mean(delta_Cb(:,:,:,:),3),4),2));

      theta_Cb = squeeze(ersp_data1_erp(:,freqs>=6&freqs<=10,:,:));
      theta_Cb = squeeze(mean(mean(mean(theta_Cb(:,:,:,:),3),4),2));

    save([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_Cb.mat'],'ersp_data1_erp','ersp_data1_cmr',...
          'delta_Cb','theta_Cb');
%       if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_M1.mat'],'file')
%             load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_M1.mat']);
%       end
    clear ersp_data1_erp ersp_data1_cmr delta_Cb theta_Cb data

      end
    end
end
toc;
%% Read collated LFP file for M1 and analyze

clear;clc;close all; tic;

bmi_session_info;

for s=M1
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Collated_LFP_new1.mat'],'file') && exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_M1.mat'],'file')

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Collated_LFP_new1.mat'],'badchans1');
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trials');

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_M1.mat']);

        
       % M1
        data = lfp_M1;
        data(badchans1,:)  = [];
        clear lfp_M1
        event = all_trials(1);
        % get 2 min spontaneous activity before session start
        if(event-60*Fs > 0)
            data = data(:,event-60*Fs:event);

      
      % Z-scoring
      for j=1:size(data,1)
        tmp = reshape(squeeze(data(j,:)),1,[]);
        trial_data_n(j,:,:) = (data(j,:,:)-mean(tmp))./std(tmp);
      end
      
      % Median subtraction
      med = median(trial_data_n,1);
      trial_data1_cmr = bsxfun(@minus,trial_data_n,med);
      
      % Get Event=Related Spectral Perturbation (ERSP) using eeglab for M1 without ERP subtraction
      for j=1:size(trial_data1_cmr,1)
        disp(['Ch_',num2str(j)]);
        data = squeeze(trial_data1_cmr(j,:,:));
        [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,size(data,1),[0 size(data,1)-1],Fs,[0.01 0.1],...
          'baseline',NaN,'freqs', [0 10],'verbose','off','trialbase','on','plotersp','off','plotitc','off');
        clf;
      end
      
      % Get normalized ESRP data generated before subtracting ERP
      data     = ersp_data.*conj(ersp_data);
      bl       = mean(data,3);
      bl_std   = std(data,[],3);
      datanorm = bsxfun(@minus,data,bl);
      ersp_datanorm1_cmr = bsxfun(@rdivide,datanorm,bl_std);
      ersp_data1_cmr = data;
      
      % ERP subtraction (Channel by Channel)
      for ch=1:size(trial_data1_cmr,1)
        single_channel_data = squeeze(trial_data1_cmr(ch,:,:));
        mean_erp = mean(single_channel_data,1);
        trial_data1_erp(ch,:,:) = bsxfun(@minus,single_channel_data,mean_erp);
      end
      
      % Get Event=Related Spectral Perturbation (ERSP) using eeglab for M1 with ERP subtraction
      for j=1:size(trial_data1_erp,1)
        disp(['Ch_',num2str(j)]);
        data = squeeze(trial_data1_erp(j,:,:));
        [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,length(data),[1 length(data)],Fs,[0.01 0.1],...
          'baseline',NaN,'freqs', [0 10],'verbose','off','trialbase','on','plotersp','off','plotitc','off');
        clf;
      end
      %,'plotersp','off','plotitc','off'
      % Get normalized ESRP data generated after subtracting ERP
      data = ersp_data.*conj(ersp_data);
      bl       = mean(data,3);
      bl_std   = std(data,[],3);
      datanorm = bsxfun(@minus,data,bl);
      ersp_datanorm1_erp = bsxfun(@rdivide,datanorm,bl_std);
      ersp_data1_erp = data;

      delta_M1 = squeeze(ersp_data1_erp(:,freqs>=0.4&freqs<=4,:,:));
      delta_M1 = squeeze(mean(mean(mean(delta_M1(:,:,:,:),3),4),2));

      theta_M1 = squeeze(ersp_data1_erp(:,freqs>=6&freqs<=10,:,:));
      theta_M1 = squeeze(mean(mean(mean(theta_M1(:,:,:,:),3),4),2));

    save([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'ersp_data1_erp','ersp_data1_cmr',...
          'delta_M1','theta_M1');
%       if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_M1.mat'],'file')
%             load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_M1.mat']);
%       end
    clear ersp_data1_erp ersp_data1_cmr delta_M1 theta_M1 data*
    
        % Clear variables 
    clear trial_data_n ersp_data trial_data1_erp trial_data1_cmr
        end
      end
    end
end
toc;

%% Analyze delta and theta power in Cb session

clear;clc;close all; tic;

bmi_session_info;

delta = {};
theta = {};

for s=intact
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_Cb.mat'],'file')

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_Cb.mat'],'delta_Cb','theta_Cb');
            load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      
      end
      delta_intact{s,n} =  delta_Cb;
      theta_intact{s,n} = theta_Cb;
            time_intact{s,n} = mean_perf(1);
      suc_intact{s,n} = mean_perf(2);
      rat_id_intact{s,n} = s;
    end
    
end


for s=stroke
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_Cb.mat'],'file')

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_Cb.mat'],'delta_Cb','theta_Cb');
            load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      end
      delta_stroke{s,n} =  delta_Cb;
      theta_stroke{s,n} = theta_Cb;
            time_stroke{s,n} = mean_perf(1);
      suc_stroke{s,n} = mean_perf(2);
            rat_id_stroke{s,n} = s;
    end
    
end
%%
% Compare delta psd intact to stroke
delta_intact = cellfun(@mean,delta_intact);
delta_intact = delta_intact(~isnan(delta_intact));

delta_stroke = cellfun(@mean,delta_stroke);
delta_stroke = delta_stroke(~isnan(delta_stroke));
mean(delta_intact) 
mean(delta_stroke)
[h,p] = ttest2(delta_intact,delta_stroke)

bar([mean(delta_intact),mean(delta_stroke)],0.6,'w'); hold on;
errorbar([1,2],[mean(delta_intact),mean(delta_stroke)],[std(delta_intact)/sqrt((size(delta_intact,1)-1)),std(delta_stroke)/sqrt((size(delta_stroke,1)-1))],'color','black','linewidth',4);
set(gca,'TickDir','out');
% ylim([0,15]);

% Compare theta psd intact to stroke
theta_intact = cellfun(@mean,theta_intact);
theta_intact = theta_intact(~isnan(theta_intact));

theta_stroke = cellfun(@mean,theta_stroke);
theta_stroke = theta_stroke(~isnan(theta_stroke));
mean(theta_intact) 
mean(theta_stroke)
[h,p] = ttest2(theta_intact,theta_stroke)
figure;
bar([mean(theta_intact),mean(theta_stroke)],0.6,'w'); hold on;
errorbar([1,2],[mean(theta_intact),mean(theta_stroke)],[std(theta_intact)/sqrt((size(theta_intact,1)-1)),std(theta_stroke)/sqrt((size(theta_stroke,1)-1))],'color','black','linewidth',4);
set(gca,'TickDir','out');
ylim([0,20]);

%% Analyze delta and theta power in M1 session

clear;clc;close all; tic;

bmi_session_info;

delta = {};
theta = {};

for s=M1_intact
%     s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'file')
         s
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'delta_M1','theta_M1');
      end
      delta_intact{s,n} =  delta_M1;
      theta_intact{s,n} = theta_M1;
    end
    
end


for s=M1_stroke
%     s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'file')
        s
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'delta_M1','theta_M1');
      end
      delta_stroke{s,n} =  delta_M1;
      theta_stroke{s,n} = theta_M1;
    end
    
end
%
% Compare delta psd intact to stroke
delta_intact = cellfun(@mean,delta_intact);
delta_intact = delta_intact(~isnan(delta_intact));

delta_stroke = cellfun(@mean,delta_stroke);
delta_stroke = delta_stroke(~isnan(delta_stroke));
mean(delta_intact) 
mean(delta_stroke)
[h,p] = ttest2(delta_intact,delta_stroke,'Tail','left')

bar([mean(delta_intact),mean(delta_stroke)],0.6,'w'); hold on;
errorbar([1,2],[mean(delta_intact),mean(delta_stroke)],[std(delta_intact)/sqrt((size(delta_intact,1)-1)),std(delta_stroke)/sqrt((size(delta_stroke,1)-1))],'color','black','linewidth',4);
set(gca,'TickDir','out');
ylim([0,15]);
% Compare theta psd intact to stroke
theta_intact = cellfun(@mean,theta_intact);
theta_intact = theta_intact(~isnan(theta_intact));

theta_stroke = cellfun(@mean,theta_stroke);
theta_stroke = theta_stroke(~isnan(theta_stroke));
mean(theta_intact) 
mean(theta_stroke)
figure;
[h,p] = ttest2(theta_intact,theta_stroke,'Tail','left')
bar([mean(theta_intact),mean(theta_stroke)],0.6,'w'); hold on;
errorbar([1,2],[mean(theta_intact),mean(theta_stroke)],[std(theta_intact)/sqrt((size(theta_intact,1)-1)),std(theta_stroke)/sqrt((size(theta_stroke,1)-1))],'color','black','linewidth',4);
set(gca,'TickDir','out');
ylim([0,15]);


%% Analyze delta and theta power in M1 session and BMI performance correlation 

clear;clc;close all; tic;

bmi_session_info;

delta = {};
theta = {};

for s=M1_intact
%     s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'file')
         s
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'delta_M1','theta_M1');
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      end
      delta_intact{s,n} =  delta_M1;
      theta_intact{s,n} = theta_M1;
      time_intact{s,n} = mean_perf(1);
      suc_intact{s,n} = mean_perf(2);
      rat_id_intact{s,n} = s;
    end
    
end


for s=M1_stroke
%     s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1
      if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'file')
        s
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\LFP_power_spont_M1.mat'],'delta_M1','theta_M1');
            load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      end
      delta_stroke{s,n} =  delta_M1;
      theta_stroke{s,n} = theta_M1;
      time_stroke{s,n} = mean_perf(1);
      suc_stroke{s,n} = mean_perf(2);
            rat_id_stroke{s,n} = s;
    end
    
end
%%
% Compare delta psd intact to stroke
delta_intact = cellfun(@mean,delta_intact);
delta_intact = delta_intact(~isnan(delta_intact));

delta_stroke = cellfun(@mean,delta_stroke);
delta_stroke = delta_stroke(~isnan(delta_stroke));

% Compare theta psd intact to stroke
theta_intact = cellfun(@mean,theta_intact);
theta_intact = theta_intact(~isnan(theta_intact));

theta_stroke = cellfun(@mean,theta_stroke);
theta_stroke = theta_stroke(~isnan(theta_stroke));

% time and success rate for BMI 
time_intact = time_intact(~cellfun(@isempty,time_intact));
time_intact = cell2mat(time_intact);
suc_intact = suc_intact(~cellfun(@isempty,suc_intact));
suc_intact = cell2mat(suc_intact);

time_stroke = time_stroke(~cellfun(@isempty,time_stroke));
time_stroke = cell2mat(time_stroke);
suc_stroke = suc_stroke(~cellfun(@isempty,suc_stroke));
suc_stroke = cell2mat(suc_stroke);
rat_id_intact = rat_id_intact(~cellfun(@isempty,rat_id_intact));
rat_id_intact = cell2mat(rat_id_intact);
rat_id_stroke = rat_id_stroke(~cellfun(@isempty,rat_id_stroke));
rat_id_stroke = cell2mat(rat_id_stroke);
%% correlation between delta power and BMI perfromance

[r,p] = corr([delta_intact;delta_stroke],[time_intact;time_stroke])

[r,p] = corr([delta_intact;delta_stroke],[suc_intact;suc_stroke])

%%
[r,p] = corr([theta_intact;theta_stroke],[time_intact;time_stroke])

[r,p] = corr([theta_intact;theta_stroke],[suc_intact;suc_stroke])

%% Linear mixed effect - delta power


tbl = [[delta_intact;delta_stroke],[repelem(1,length(delta_intact))';repelem(2,length(delta_stroke))'],...
    [time_intact;time_stroke],[suc_intact;suc_stroke],[rat_id_intact;rat_id_stroke]];
tbl = array2table(tbl,'VariableNames',{'delta','stroke','time','success','RatID'});
formula = 'delta ~ time + success + stroke + (time + success + stroke | RatID)';
lme_R2 = fitlme(tbl,formula)

%% Linear mixed effect - theta power
tbl = [[theta_intact;theta_stroke],[repelem(1,length(theta_intact))';repelem(2,length(theta_stroke))'],...
    [time_intact;time_stroke],[suc_intact;suc_stroke],[rat_id_intact;rat_id_stroke]];
tbl = array2table(tbl,'VariableNames',{'theta','stroke','time','success','RatID'});
formula = 'theta ~ time + success + stroke + (time + success + stroke | RatID)';
lme_R2 = fitlme(tbl,formula)