%% Get Event-Related Spectral Perturbation (ERSP)
%  Author - Aamir Abbasi
% Modified - Rohit
%% -------------------------------------------------------
clear;clc;
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\I128\Data\';
% addpath(genpath('Z:\Matlab for analysis\eeglab'));

bmiBlocks = dir(rootpath);
         

for i=3:6%:9%3:length(bmiBlocks)
  
  disp(bmiBlocks(i).name);
  
  % Read collated LFP matrix
  load([rootpath, bmiBlocks(i).name,'\','Collated_LFP_new1.mat']);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR M1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Clear variables common with M1
%   clear trial_data_n ersp_data trial_data1_erp 
%   
%   % Z-scoring
%   for j=1:size(trial_data1,1)
%     tmp = reshape(squeeze(trial_data1(j,:,:)),1,[]);
%     trial_data_n(j,:,:) = (trial_data1(j,:,:)-mean(tmp))./std(tmp);
%   end
%   
%   % Median subtraction
%   med = median(trial_data_n,1);
%   trial_data1_cmr = bsxfun(@minus,trial_data_n,med);
%   
%   % Get Event=Related Spectral Perturbation (ERSP) using eeglab for M1 without ERP subtraction
%   for j=1:size(trial_data1_cmr,1)
%     disp(['Ch_',num2str(j)]);
%     data = squeeze(trial_data1_cmr(j,1:8000,:));
%     [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,8000,[-4000 4000],Fs,[0.01 0.1],...
%       'baseline',NaN,'freqs', [0 60],'verbose','off','trialbase','on');
%     clf;
%   end
%   
%   % Get normalized ESRP data generated before subtracting ERP
%   data     = ersp_data.*conj(ersp_data);
%   bl       = mean(data(:,:,times>-1500&times<-500,:),3);
%   bl_std   = std(data(:,:,times>-1500&times<-500,:),[],3);
%   datanorm = bsxfun(@minus,data,bl);
%   ersp_datanorm1_cmr = bsxfun(@rdivide,datanorm,bl_std);
%   ersp_data1_cmr = data;
%   
%   % ERP subtraction (Channel by Channel)
%   for ch=1:size(trial_data1_cmr,1)
%     single_channel_data = squeeze(trial_data1_cmr(ch,:,:));
%     mean_erp = mean(single_channel_data,2);
%     trial_data1_erp(ch,:,:) = bsxfun(@minus,single_channel_data,mean_erp);
%   end
%   
%   % Get Event=Related Spectral Perturbation (ERSP) using eeglab for M1 with ERP subtraction
%   for j=1:size(trial_data1_erp,1)
%     disp(['Ch_',num2str(j)]);
%     data = squeeze(trial_data1_erp(j,1:8000,:));
%     [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,8000,[-4000 4000],Fs,[0.01 0.1],...
%       'baseline',NaN,'freqs', [0 60],'verbose','off','trialbase','on');
%     clf;
%   end
%   
%   % Get normalized ESRP data generated after subtracting ERP
%   data = ersp_data.*conj(ersp_data);
%   bl       = mean(data(:,:,times>-1500&times<-500,:),3);
%   bl_std   = std(data(:,:,times>-1500&times<-500,:),[],3);
%   datanorm = bsxfun(@minus,data,bl);
%   ersp_datanorm1_erp = bsxfun(@rdivide,datanorm,bl_std);
%   ersp_data1_erp = data;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR CEREBELLUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Clear variables common with M1
  clear trial_data_n ersp_data trial_data2_erp
  
  % Z-scoring
  for j=1:size(trial_data2,1)
    tmp = reshape(squeeze(trial_data2(j,:,:)),1,[]);
    trial_data_n(j,:,:) = (trial_data2(j,:,:)-mean(tmp))./std(tmp);
  end
  
  % Median subtraction CMR
  med = median(trial_data_n,1);
  trial_data2_cmr = bsxfun(@minus,trial_data_n,med);
  
  % Get Event=Related Spectral Perturbation (ERSP) using eeglab for Cb without ERP subtraction
  for j=1:size(trial_data2_cmr,1)
    disp(['Ch_',num2str(j)]);
    data = squeeze(trial_data2_cmr(j,1:8000,:));
    [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,8000,[-4000 4000],Fs,[0.01 0.1],...
      'baseline',NaN,'freqs', [0 60],'verbose','off','trialbase','on');
    clf;
  end
  
  % Get normalized ESRP data generated before subtracting ERP
  data = ersp_data.*conj(ersp_data);
  bl=mean(data(:,:,times>-1500&times<-500,:),3);
  bl_std=std(data(:,:,times>-1500&times<-500,:),[],3);
  datanorm=bsxfun(@minus,data,bl);
  ersp_datanorm2_cmr = bsxfun(@rdivide,datanorm,bl_std);
  ersp_data2_cmr = data;
  
  % ERP subtraction (Channel by Channel)
  for ch=1:size(trial_data2_cmr,1)
    single_channel_data = squeeze(trial_data2_cmr(ch,:,:));
    mean_erp = mean(single_channel_data,2);
    trial_data2_erp(ch,:,:) = bsxfun(@minus,single_channel_data,mean_erp);
  end
  
  % Get Event=Related Spectral Perturbation (ERSP) using eeglab for Cb with ERP subtraction
  for j=1:size(trial_data2_erp,1)
    disp(['Ch_',num2str(j)]);
    data = squeeze(trial_data2_erp(j,1:8000,:));
    [~,~,~,times,freqs,~,~,ersp_data(j,:,:,:)] = newtimef(data,8000,[-4000 4000],Fs,[0.01 0.1],...
      'baseline',NaN,'freqs', [0 60],'verbose','off','trialbase','on');
    clf;
  end
  
  % Get normalized ESRP data generated after subtracting ERP
  data = ersp_data.*conj(ersp_data);
  bl=mean(data(:,:,times>-1500&times<-500,:),3);
  bl_std=std(data(:,:,times>-1500&times<-500,:),[],3);
  datanorm=bsxfun(@minus,data,bl);
  ersp_datanorm2_erp = bsxfun(@rdivide,datanorm,bl_std);
  ersp_data2_erp = data;
  
  % Save
  save([rootpath, bmiBlocks(i).name,'\','ERSP_Norm.mat'],...%'ersp_datanorm1_cmr','ersp_datanorm1_erp',...
    'ersp_datanorm2_cmr','ersp_datanorm2_erp','times','freqs','-v7.3');  
  save([rootpath, bmiBlocks(i).name,'\','ERSP.mat'],...%'ersp_data1_cmr','ersp_data1_erp',...
    'ersp_data2_cmr','ersp_data2_erp','times','freqs','-v7.3');   
end
close;
disp('done');