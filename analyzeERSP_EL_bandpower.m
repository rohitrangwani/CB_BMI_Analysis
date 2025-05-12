%% THIS SCRIPT IS USED FOR POWER ANALYSIS OF EARLY AND LATE TRIALS
%- Author @Aamir Abbasi
%-------------------------------------------------------------------------------------------
%% Get all session data for early and late trials
%  ----------------------------------------------------------------------------------------------------------------------
clear; clc; close all;
disp('running...');
% rootpath = 'Z:\Rohit\BMI_Data\I122\Data\';
   

    meT_M1 = [];
    mlT_M1 = [];
    meT_Cb = [];
    mlT_Cb = [];

freq_band = {[0.1 4] [3 6] [6 14] [13 30]};

sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};

stroke= [4:5,9:10];

intact = [1:10];
intact = setdiff(intact, stroke);

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]...
    [4,6,8],[5:9],[8,10],[7:9],[3,5,6]};
path = 'Z:\Rohit\BMI_Data\';


for s =2
    s
    
  rootpath = ['Z:\Rohit\BMI_Data\',sub{s},'\Data\'];  
  % Display current block id 
  bmiBlocks = dir(rootpath); 
  
for f =1:length(freq_band)
    
  meT_M1 = [];
  mlT_M1 = [];
  meT_Cb = [];
  mlT_Cb = [];
  for i=robust_session{s}%4:length(bmiBlocks) %
  bmiBlocks(i).name
  
  % Load raw LFP channels
  load([rootpath, bmiBlocks(i).name,'\','Collated_LFP_CMR_ERP.mat'],'trial_data1_erp','trial_data2_erp','Fs','valid_performance'); 
  
  ersp_data1_erp = [];
  % Load ERSP data for current block
  load([rootpath, bmiBlocks(i).name,'\','ERSP.mat'],'ersp_data1_erp','ersp_data2_erp','times','freqs');  

  
  % Adjust trial start time
  times = times - 50;  
  
  % Store M1 and Cb ERSP Data for BMI1
  data1 = ersp_data1_erp;
  data2 = ersp_data2_erp;

  ntrialsBMI1 = size(data1,4); 
  ntrialsBMI2 = size(data2,4);
  % Normalize ERSP
  if(data1)
  data = data1;
      bl   = mean(data(:,:,times>-1500&times<500,:),3);  
      bl_std = std(data(:,:,times>-1500&times<500,:),[],3);  
      datanorm = bsxfun(@minus,data,bl);  
      data1 = bsxfun(@rdivide,datanorm,bl_std);  
  end
  data = data2;
  bl   = mean(data(:,:,times>-1500&times<500,:),3);  
  bl_std = std(data(:,:,times>-1500&times<500,:),[],3);  
  datanorm = bsxfun(@minus,data,bl);  
  data2 = bsxfun(@rdivide,datanorm,bl_std);   

  ersp_data1_erp = data1(:,:,:,1:ntrialsBMI1);
  ersp_data2_erp = data2(:,:,:,1:ntrialsBMI2);
  
  if(data1) 
  % Get M1 early/late trials ersp
  eT_M1 = squeeze(ersp_data1_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,1:floor(size(ersp_data1_erp,4)/3)));
  lT_M1 = squeeze(ersp_data1_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,end-floor(size(ersp_data1_erp,4)/3)+1:end));
  
  % Get only slow trials from early 
  eT_M1 = squeeze(eT_M1(:,:,:,logical((valid_performance(1:floor(size(trial_data2_erp,3)/3))<=15).*...
    (valid_performance(1:floor(size(trial_data2_erp,3)/3))>=5))));

  % Get only fast trials from late
  lT_M1 = squeeze(lT_M1(:,:,:,logical((valid_performance(end-floor(size(trial_data2_erp,3)/3)+1:end)<=5))));  

  % Mean across trials/freqs/time
  eT_M1 = squeeze(mean(mean(mean(eT_M1(:,:,times>=200&times<=600,:),3),4),2));
  lT_M1 = squeeze(mean(mean(mean(lT_M1(:,:,times>=200&times<=600,:),3),4),2));
  early_M1 = eT_M1(lT_M1-eT_M1>0);
  late_M1  = lT_M1(lT_M1-eT_M1>0);
  
  % Save channels showing an increase in power
  m1_chans = lT_M1-eT_M1>0;  

  meT_M1 = [meT_M1;nanmean(early_M1)];
  mlT_M1 = [mlT_M1;nanmean(late_M1)];

  end
  % Get Cb early/late trials ersp 
  eT_Cb = squeeze(ersp_data2_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,1:floor(size(ersp_data2_erp,4)/3)));
  lT_Cb = squeeze(ersp_data2_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,end-floor(size(ersp_data2_erp,4)/3)+1:end));
  
  % Get only slow trials from early  
  eT_Cb = squeeze(eT_Cb(:,:,:,logical((valid_performance(1:floor(size(trial_data2_erp,3)/3))<=15).*...
    (valid_performance(1:floor(size(trial_data2_erp,3)/3))>=5))));  
  
  % Get only fast trials from late
  lT_Cb = squeeze(lT_Cb(:,:,:,logical((valid_performance(end-floor(size(trial_data2_erp,3)/3)+1:end)<=5))));  
  
  % Mean across trials/freqs/time
  eT_Cb = squeeze(mean(mean(mean(eT_Cb(:,:,times>=200&times<=600,:),3),4),2));
  lT_Cb = squeeze(mean(mean(mean(lT_Cb(:,:,times>=200&times<=600,:),3),4),2));  
  early_Cb = eT_Cb(lT_Cb-eT_Cb>0);
  late_Cb  = lT_Cb(lT_Cb-eT_Cb>0);
  
  % Save channels showing an increase in power
% %   m1_chans = [];
  cb_chans = lT_Cb-eT_Cb>0; 

%   save([rootpath,bmiBlocks(i).name,'\ERSP_chans.mat'],'m1_chans','cb_chans');
  
  % Store average PSD across sessions
%   meT_M1 = [meT_M1;nanmean(early_M1) str2num(bmiBlocks(i)(3:4))];
%   mlT_M1 = [mlT_M1;nanmean(late_M1) str2num(bmiBlocks(i)(3:4))];
  
%   meT_Cb = [meT_Cb;nanmean(early_Cb) str2num(bmiBlocks(i).name(3:4))];
%   mlT_Cb = [mlT_Cb;nanmean(late_Cb) str2num(bmiBlocks(i).name(3:4))];

  meT_Cb = [meT_Cb;nanmean(early_Cb)];
  mlT_Cb = [mlT_Cb;nanmean(late_Cb)];
  

  end  
  
  
  
  % Replace all negative vals with zeros
meT_M1(meT_M1<0) = 0;
mlT_M1(mlT_M1<0) = 0;
meT_Cb(meT_Cb<0) = 0;
mlT_Cb(mlT_Cb<0) = 0;

% Z-score PSD
% meT_M1 = zscore(meT_M1')';
% mlT_M1 = zscore(mlT_M1')';
% meT_Cb = zscore(meT_Cb')';
% mlT_Cb = zscore(mlT_Cb')';

    bar(nanmean([meT_M1,mlT_M1])); hold on;
%     plot([1,meT_Cb(:)],[2,mlT_Cb(:)]);
display(meT_M1)
    for j = 1:length(meT_M1)
        
        plot([[meT_M1(j)];[mlT_M1(j)]],'k');hold on;
    end
    errorbar(1,nanmean(meT_M1),nanstd(meT_M1)./sqrt(length(meT_M1)-1));
    errorbar(2,nanmean(mlT_M1),nanstd(mlT_M1)./sqrt(length(mlT_M1)-1));
    [h_M1,p_M1] = ttest(meT_M1,mlT_M1);
    display(p_M1);
    
    set(gcf,'WindowState','maximized')
    saveas(gcf,[rootpath,'\..\BandPower_M1_',num2str(f),'.tiff']);
    close;
    
    
 %     subplot(length(freq_band),1,f)
    bar(nanmean([meT_Cb,mlT_Cb])); hold on;
%     plot([1,meT_Cb(:)],[2,mlT_Cb(:)]);
    display(meT_Cb)
    for j = 1:length(meT_Cb)
        
        plot([[meT_Cb(j)];[mlT_Cb(j)]],'k');hold on;
    end
    errorbar(1,nanmean(meT_Cb),nanstd(meT_Cb)./sqrt(length(meT_Cb)-1));
    errorbar(2,nanmean(mlT_Cb),nanstd(mlT_Cb)./sqrt(length(mlT_Cb)-1));
    [h_Cb,p_Cb] = ttest(meT_Cb,mlT_Cb);
    display(p_Cb);

    save([rootpath,'\..\LFP_power_consolidated',num2str(f),'.mat'],'meT_M1','mlT_M1',...
          'meT_Cb','mlT_Cb');
    
    set(gcf,'WindowState','maximized')
    saveas(gcf,[rootpath,'\..\BandPower_Cb_',num2str(f),'.tiff']);
    close;
end
end


% % Replace all negative vals with zeros
% meT_M1(meT_M1<0) = 0;
% mlT_M1(mlT_M1<0) = 0;
% meT_Cb(meT_Cb<0) = 0;
% mlT_Cb(mlT_Cb<0) = 0;
% 
% % Z-score PSD
% meT_M1 = zscore(meT_M1')';
% mlT_M1 = zscore(mlT_M1')';
% meT_Cb = zscore(meT_Cb')';
% mlT_Cb = zscore(mlT_Cb')';
% save([rootpath,'\LFP_power_consolidated.mat'],'meT_M1','mlT_M1',...
%       'meT_Cb','mlT_Cb');
  
% save([rootpath,'\LFP_power_consolidated.mat'],'meT_Cb','mlT_Cb');
% disp('done!');


%% Get all session data for early and late trials - Coherence analysis - FFC 
%  ----------------------------------------------------------------------------------------------------------------------
clear; clc; close all;
disp('running...');

freq_band = {[0.1 4] [3 6] [6 14] [13 30]};
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};

stroke= [4:5,9:10];

intact = [1:10];
intact = setdiff(intact, stroke);

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]...
    [4,6,8],[5:9],[10],[7:9],[3,5,6]};
path = 'Z:\Rohit\BMI_Data\';

M1 = [1,2,6,7,8];

for s =8%M1
    s
  rootpath = ['Z:\Rohit\BMI_Data\',sub{s},'\Data\'];  
  % Display current block id 
  bmiBlocks = dir(rootpath); 
for f =4%1:length(freq_band)
    
  meT = [];
  mlT = [];
  for i=robust_session{s}%length(bmiBlocks) %
    
  % Display current block id 
  disp(bmiBlocks(i).name);
  
  % Load raw LFP channels
%   load([rootpath, bmiBlocks(i).name,'\','Collated_LFP_CMR_ERP.mat'],'trial_data1_erp','trial_data2_erp','Fs','valid_performance'); 
  
%   % Load ERSP data for current block
  load([rootpath, bmiBlocks(i).name,'\','ERSP.mat'],'ersp_data1_erp','ersp_data2_erp','times','freqs');  
  
  % load coherence data
  load([rootpath,bmiBlocks(i).name,'\FFC_EL_ERP.mat'],'e_coher_erp','l_coher_erp' );

  
  % Adjust trial start time
%   times = times - 50;  
  
  % Store M1 and Cb ERSP Data for BMI1
  data1 = e_coher_erp;
  data2 = l_coher_erp;
%   ntrialsBMI1 = size(data2,4);
   
  % Normalize ERSP
%   data = data1;
%   bl   = mean(data(:,:,times>-3000&times<-500,:),3);  
%   bl_std = std(data(:,:,times>-3000&times<-500,:),[],3);  
%   datanorm = bsxfun(@minus,data,bl);  
%   data1 = bsxfun(@rdivide,datanorm,bl_std);  
%   
%   data = data2;
%   bl   = mean(data(:,:,times>-3000&times<-500,:),3);  
%   bl_std = std(data(:,:,times>-3000&times<-500,:),[],3);  
%   datanorm = bsxfun(@minus,data,bl);  
%   data2 = bsxfun(@rdivide,datanorm,bl_std);   

%   ersp_data1_erp = data1(:,:,:,1:ntrialsBMI1);
%   ersp_data2_erp = data2(:,:,:,1:ntrialsBMI1);
  
  
  % Get M1 early/late trials ersp
  eT = squeeze(e_coher_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,:));
  lT = squeeze(l_coher_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,:));

  
  % Mean across trials/freqs/time
  eT = squeeze(mean(mean(mean(eT(:,:,:,:),3),4),2));
  lT = squeeze(mean(mean(mean(lT(:,:,:,:),3),4),2));
  early = eT(lT-eT>0);
  late  = lT(lT-eT>0);  
  
  
%   % Save channels showing an increase in power
%   m1_chans = lT_M1-eT_M1>0;
% % %   m1_chans = [];
%   cb_chans = lT_Cb-eT_Cb>0; 
% 
%   save([rootpath,bmiBlocks(i).name,'\ERSP_chans.mat'],'m1_chans','cb_chans');
  
  % Store average PSD across sessions
%   meT_M1 = [meT_M1;nanmean(early_M1) str2num(bmiBlocks(i)(3:4))];
%   mlT_M1 = [mlT_M1;nanmean(late_M1) str2num(bmiBlocks(i)(3:4))];
  
%   meT_Cb = [meT_Cb;nanmean(early_Cb) str2num(bmiBlocks(i).name(3:4))];
%   mlT_Cb = [mlT_Cb;nanmean(late_Cb) str2num(bmiBlocks(i).name(3:4))];
  meT = [meT;nanmean(early)];
  mlT = [mlT;nanmean(late)];

  end  

  % Replace all negative vals with zeros
meT(meT<0) = 0;
mlT(mlT<0) = 0;


% Z-score PSD
% meT_M1 = zscore(meT_M1')';
% mlT_M1 = zscore(mlT_M1')';
% meT_Cb = zscore(meT_Cb')';
% mlT_Cb = zscore(mlT_Cb')';

    bar(nanmean([meT,mlT])); hold on;
%     plot([1,meT_Cb(:)],[2,mlT_Cb(:)]);
display(meT)
    for j = 1:length(meT)
        
        plot([[meT(j)];[mlT(j)]],'k');hold on;
    end
    errorbar(1,nanmean(meT),nanstd(meT)./sqrt(length(meT)-1));
    errorbar(2,nanmean(mlT),nanstd(mlT)./sqrt(length(mlT)-1));
    [h_M1,p_M1] = ttest(meT,mlT);
    display(p_M1);
    
    set(gcf,'WindowState','maximized')
    saveas(gcf,[rootpath,'\..\BandPower_FFC_',num2str(f),'.tiff']);
    close;
    
   
    save([rootpath,'\..\FFC_power_consolidated',num2str(f),'.mat'],'meT','mlT');
    
end
end    
% % Replace all negative vals with zeros
% meT_M1(meT_M1<0) = 0;
% mlT_M1(mlT_M1<0) = 0;
% meT_Cb(meT_Cb<0) = 0;
% mlT_Cb(mlT_Cb<0) = 0;
% 
% % Z-score PSD
% meT_M1 = zscore(meT_M1')';
% mlT_M1 = zscore(mlT_M1')';
% meT_Cb = zscore(meT_Cb')';
% mlT_Cb = zscore(mlT_Cb')';
% save([rootpath,'\LFP_power_consolidated.mat'],'meT_M1','mlT_M1',...
%       'meT_Cb','mlT_Cb');
  
% save([rootpath,'\LFP_power_consolidated.mat'],'meT_Cb','mlT_Cb');
% disp('done!');


%% Get all session data for early and late trials - Cb only
%  ----------------------------------------------------------------------------------------------------------------------
clear; clc; close all;
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\I112\Data\';
bmiBlocks = dir(rootpath);
 
freq_band = {[0.1 4] [3 6] [6 14]};

for f =1:length(freq_band)
  meT_M1 = [];
  mlT_M1 = [];    
  meT_Cb = [];
  mlT_Cb = [];
  for i=[5:7,9:10]%[3:5 8:length(bmiBlocks)]%[5 7:10 13:14]%3:6%length(bmiBlocks) %
    
  % Display current block id 
  disp(bmiBlocks(i).name);
  
  % Load raw LFP channels
  load([rootpath, bmiBlocks(i).name,'\','Collated_LFP_CMR_ERP.mat'],'trial_data1_erp','trial_data2_erp','Fs','valid_performance'); 
  
  % Load ERSP data for current block
  load([rootpath, bmiBlocks(i).name,'\','ERSP.mat'],'ersp_data1_erp','ersp_data2_erp','times','freqs');  

  
  % Adjust trial start time
  times = times - 50;  
  
  % Store M1 and Cb ERSP Data for BMI1
  data2 = ersp_data2_erp;
  ntrialsBMI1 = size(data2,4);
   
  % Normalize ERSP
%   data = data1;
%   bl   = mean(data(:,:,times>-3000&times<-500,:),3);  
%   bl_std = std(data(:,:,times>-3000&times<-500,:),[],3);  
%   datanorm = bsxfun(@minus,data,bl);  
%   data1 = bsxfun(@rdivide,datanorm,bl_std);  
  
  data = data2;
  bl   = mean(data(:,:,times>-3000&times<-500,:),3);  
  bl_std = std(data(:,:,times>-3000&times<-500,:),[],3);  
  datanorm = bsxfun(@minus,data,bl);  
  data2 = bsxfun(@rdivide,datanorm,bl_std);   

  ersp_data2_erp = data2(:,:,:,1:ntrialsBMI1);

  % Get Cb early/late trials ersp 
  eT_Cb = squeeze(ersp_data2_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,1:floor(size(ersp_data2_erp,4)/3)));
  lT_Cb = squeeze(ersp_data2_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,end-floor(size(ersp_data2_erp,4)/3)+1:end));
  
   % Get only slow trials from early 
  eT_Cb = squeeze(eT_Cb(:,:,:,logical((valid_performance(1:floor(size(trial_data2_erp,3)/3))<=15).*...
    (valid_performance(1:floor(size(trial_data2_erp,3)/3))>=5))));  
  % Get only fast trials from late 
  lT_Cb = squeeze(lT_Cb(:,:,:,logical((valid_performance(end-floor(size(trial_data2_erp,3)/3)+1:end)<=5))));  

  
  eT_Cb = squeeze(mean(mean(mean(eT_Cb(:,:,times>=200&times<=600,:),3),4),2));
  lT_Cb = squeeze(mean(mean(mean(lT_Cb(:,:,times>=200&times<=600,:),3),4),2));  
  early_Cb = eT_Cb(lT_Cb-eT_Cb>0);
  late_Cb  = lT_Cb(lT_Cb-eT_Cb>0);
  
  % Save channels showing an increase in power
  m1_chans = [];
  cb_chans = lT_Cb-eT_Cb>0; 
  save([rootpath,bmiBlocks(i).name,'\ERSP_chans.mat'],'m1_chans','cb_chans');
  
  % Store average PSD across sessions
  meT_Cb = [meT_Cb;nanmean(early_Cb)];
  mlT_Cb = [mlT_Cb;nanmean(late_Cb)];
  

  end  
 
  
    % Replace all negative vals with zeros
meT_M1(meT_M1<0) = 0;
mlT_M1(mlT_M1<0) = 0;
meT_Cb(meT_Cb<0) = 0;
mlT_Cb(mlT_Cb<0) = 0;

% Z-score PSD
% meT_M1 = zscore(meT_M1')';
% mlT_M1 = zscore(mlT_M1')';
% meT_Cb = zscore(meT_Cb')';
% mlT_Cb = zscore(mlT_Cb')';

 %     subplot(length(freq_band),1,f)
    bar(nanmean([meT_Cb,mlT_Cb])); hold on;
%     plot([1,meT_Cb(:)],[2,mlT_Cb(:)]);
    display(meT_Cb)
    for j = 1:length(meT_Cb)
        
        plot([[meT_Cb(j)];[mlT_Cb(j)]],'k');hold on;
    end
    errorbar(1,nanmean(meT_Cb),nanstd(meT_Cb)./sqrt(length(meT_Cb)-1));
    errorbar(2,nanmean(mlT_Cb),nanstd(mlT_Cb)./sqrt(length(mlT_Cb)-1));
    [h_Cb,p_Cb] = ttest(meT_Cb,mlT_Cb);
    display(p_Cb);

    save([rootpath,'\LFP_power_consolidated',num2str(f),'.mat'],'meT_M1','mlT_M1',...
          'meT_Cb','mlT_Cb');
    
    set(gcf,'WindowState','maximized')
    saveas(gcf,[rootpath,'BandPower_Cb_',num2str(f),'.tiff']);
    close;
end
    
% % Replace all negative vals with zeros
% meT_M1(meT_M1<0) = 0;
% mlT_M1(mlT_M1<0) = 0;
% meT_Cb(meT_Cb<0) = 0;
% mlT_Cb(mlT_Cb<0) = 0;
% 
% % Z-score PSD
% meT_M1 = zscore(meT_M1')';
% mlT_M1 = zscore(mlT_M1')';
% meT_Cb = zscore(meT_Cb')';
% mlT_Cb = zscore(mlT_Cb')';
% save([rootpath,'\LFP_power_consolidated.mat'],'meT_M1','mlT_M1',...
%       'meT_Cb','mlT_Cb');
  
% save([rootpath,'\LFP_power_consolidated.mat'],'meT_Cb','mlT_Cb');
% disp('done!');

%% Get all subject data for early and late trials  - Cb only
%  ----------------------------------------------------------------------------------------------------------------------
clear; clc; close all;
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';


sub = {'I096','I107','I110','I111','I112'};
disp('running...');

first = [4,3,4,4,3];
last = [6,6,10,10,10];

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]};

for s=1:3%4:5
bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
 
freq_band = {[0.1 4] [3 6] [6 14]};

for f =2%1:length(freq_band)
  meT_M1 = [];
  mlT_M1 = [];    
  meT_Cb = [];
  mlT_Cb = [];
  for i=robust_session{s}
    
  % Display current block id 
  disp(bmiBlocks(i).name);
  
  % Load raw LFP channels
  load([rootpath, sub{s},'\Data\', bmiBlocks(i).name,'\Collated_LFP_CMR_ERP.mat'],'trial_data1_erp','trial_data2_erp','Fs','valid_performance'); 
  
  % Load ERSP data for current block
  load([rootpath, sub{s},'\Data\', bmiBlocks(i).name,'\ERSP.mat'],'ersp_data1_erp','ersp_data2_erp','times','freqs');  

  
  % Adjust trial start time
  times = times - 50;  
  
  % Store M1 and Cb ERSP Data for BMI1
  data2 = ersp_data2_erp;
  ntrialsBMI1 = size(data2,4);
   
  % Normalize ERSP
%   data = data1;
%   bl   = mean(data(:,:,times>-3000&times<-500,:),3);  
%   bl_std = std(data(:,:,times>-3000&times<-500,:),[],3);  
%   datanorm = bsxfun(@minus,data,bl);  
%   data1 = bsxfun(@rdivide,datanorm,bl_std);  
  
  data = data2;
  bl   = mean(data(:,:,times>-3000&times<-500,:),3);  
  bl_std = std(data(:,:,times>-3000&times<-500,:),[],3);  
  datanorm = bsxfun(@minus,data,bl);  
  data2 = bsxfun(@rdivide,datanorm,bl_std);   

  ersp_data2_erp = data2(:,:,:,1:ntrialsBMI1);

  % Get Cb early/late trials ersp 
  eT_Cb = squeeze(ersp_data2_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,1:floor(size(ersp_data2_erp,4)/3)));
  lT_Cb = squeeze(ersp_data2_erp(:,freqs>=freq_band{f}(1)&freqs<=freq_band{f}(2),:,end-floor(size(ersp_data2_erp,4)/3)+1:end));
  
   % Get only slow trials from early 
  eT_Cb = squeeze(eT_Cb(:,:,:,logical((valid_performance(1:floor(size(trial_data2_erp,3)/3))<=15).*...
    (valid_performance(1:floor(size(trial_data2_erp,3)/3))>=5))));  
  % Get only fast trials from late 
  lT_Cb = squeeze(lT_Cb(:,:,:,logical((valid_performance(end-floor(size(trial_data2_erp,3)/3)+1:end)<=5))));  

  
  eT_Cb = squeeze(mean(mean(mean(eT_Cb(:,:,times>=200&times<=600,:),3),4),2));
  lT_Cb = squeeze(mean(mean(mean(lT_Cb(:,:,times>=200&times<=600,:),3),4),2));  
  early_Cb = eT_Cb(lT_Cb-eT_Cb>0);
  late_Cb  = lT_Cb(lT_Cb-eT_Cb>0);
  
  % Save channels showing an increase in power
  m1_chans = [];
  cb_chans = lT_Cb-eT_Cb>0; 
%   save([rootpath,bmiBlocks(i).name,'\ERSP_chans.mat'],'m1_chans','cb_chans');
  
  % Store average PSD across sessions
  meT_Cb = [meT_Cb;nanmean(early_Cb)];
  mlT_Cb = [mlT_Cb;nanmean(late_Cb)];
  
  meT_Cb_all{s} = meT_Cb;
  mlT_Cb_all{s} = mlT_Cb; 
  end  
 
  
    % Replace all negative vals with zeros
% meT_M1(meT_M1<0) = 0;
% mlT_M1(mlT_M1<0) = 0;
% meT_Cb(meT_Cb<0) = 0;
% mlT_Cb(mlT_Cb<0) = 0;

% Z-score PSD
% meT_M1 = zscore(meT_M1')';
% mlT_M1 = zscore(mlT_M1')';
% meT_Cb = zscore(meT_Cb')';
% mlT_Cb = zscore(mlT_Cb')';

 %     subplot(length(freq_band),1,f)
%     bar(nanmean([meT_Cb,mlT_Cb])); hold on;
% %     plot([1,meT_Cb(:)],[2,mlT_Cb(:)]);
%     display(meT_Cb)
%     for j = 1:length(meT_Cb)
%         
%         plot([[meT_Cb(j)];[mlT_Cb(j)]],'k');hold on;
%     end
%     errorbar(1,nanmean(meT_Cb),nanstd(meT_Cb)./sqrt(length(meT_Cb)-1));
%     errorbar(2,nanmean(mlT_Cb),nanstd(mlT_Cb)./sqrt(length(mlT_Cb)-1));
%     [h_Cb,p_Cb] = ttest(meT_Cb,mlT_Cb);
%     display(p_Cb);
% 
%     save([rootpath,'\LFP_power_consolidated',num2str(f),'.mat'],'meT_M1','mlT_M1',...
%           'meT_Cb','mlT_Cb');
%     
%     set(gcf,'WindowState','maximized')
%     saveas(gcf,[rootpath,'BandPower_Cb_',num2str(f),'.tiff']);
%     close;
end

end
%%

meT_Cb = cell2mat(meT_Cb_all(:));
mlT_Cb = cell2mat(mlT_Cb_all(:));
 bar(nanmean([meT_Cb,mlT_Cb])); hold on;
%     plot([1,meT_Cb(:)],[2,mlT_Cb(:)]);
    display(meT_Cb)
    for j = 1:length(meT_Cb)
        
        plot([[meT_Cb(j)];[mlT_Cb(j)]],'k');hold on;
    end
    errorbar(1,nanmean(meT_Cb),nanstd(meT_Cb)./sqrt(length(meT_Cb)-1));
    errorbar(2,nanmean(mlT_Cb),nanstd(mlT_Cb)./sqrt(length(mlT_Cb)-1));
    [h_Cb,p_Cb] = ttest(meT_Cb,mlT_Cb);
    display(p_Cb);
    
% scatter(ones(length(meT_Cb),1),meT_Cb); hold on;
% scatter(ones(length(mlT_Cb),1)*2,mlT_Cb);

%     save([rootpath,'\LFP_power_consolidated',num2str(f),'.mat'],'meT_M1','mlT_M1',...
%           'meT_Cb','mlT_Cb');
%     
%     set(gcf,'WindowState','maximized')
%     saveas(gcf,[rootpath,'BandPower_Cb_stroke','.tiff']);
%     close;
    
    
%% Bar plot
rootpath = 'Z:\Rohit\BMI_Data\I110\Data\';
load([rootpath,'\LFP_power_consolidated.mat']);

% subplot(1,2,1); bar(nanmean([meT_M1,mlT_M1])); hold on;
% errorbar(1,nanmean(meT_M1),nanstd(meT_M1)./sqrt(length(meT_M1)-1));
% errorbar(2,nanmean(mlT_M1),nanstd(mlT_M1)./sqrt(length(mlT_M1)-1));
% [h_M1,p_M1] = ttest(meT_M1,mlT_M1);
% subplot(1,2,2); 

%
meT_Cb = meT_Cb(:,1);
mlT_Cb = mlT_Cb(:,1);

bar(nanmean([meT_Cb,mlT_Cb])); hold on;
errorbar(1,nanmean(meT_Cb),nanstd(meT_Cb)./sqrt(length(meT_Cb)-1));
errorbar(2,nanmean(mlT_Cb),nanstd(mlT_Cb)./sqrt(length(mlT_Cb)-1));
[h_Cb,p_Cb] = ttest(meT_Cb,mlT_Cb);
display(p_Cb);


%% Plot session by session PSD
clc;
figure('Color','white','Position',[404 276 1439 613]);
for i=1:size(meT_Cb)
%   subplot(1,2,1); plot(freqs,meT_M1(i,:));
%   hold on; plot(freqs,mlT_M1(i,:));
%   xlim([0 15]); %set(gca,'YScale','log');
%   
%   subplot(1,2,2);
  subplot(1,2,1);plot(freqs,meT_Cb(i));
  hold on; subplot(1,2,2); plot(freqs,mlT_Cb(i));
  xlim([0 15]); %set(gca,'YScale','log');
  
  suptitle(bmiBlocks(i).name(10:end));
  pause;
%   clf
end

% Plot average PSD
% subplot(1,2,1); plot(freqs,mean(meT_M1));
% hold on; plot(freqs,mean(mlT_M1));
% xlim([0 15]); %set(gca,'YScale','log');
% 
% subplot(1,2,2);
plot(freqs,mean(meT_Cb));
hold on; plot(freqs,mean(mlT_Cb));
xlim([0 15]); %set(gca,'YScale','log');
%% Plot session by session PSD

figure('Color','white','Position',[404 276 1439 613]);
for i=1:size(meT_M1,1)
  subplot(1,2,1); plot(freqs,meT_M1(i));
  hold on; plot(freqs,mlT_M1(i,:));
  xlim([0 15]); %set(gca,'YScale','log');
  
  subplot(1,2,2); plot(freqs,meT_Cb(i));
  hold on; plot(freqs,mlT_Cb(i,:));
  xlim([0 15]); %set(gca,'YScale','log');
  
  
  suptitle(bmiBlocks(i).name(10:end));
  pause;
%   clf
end

% Plot average PSD
subplot(1,2,1); plot(freqs,mean(meT_M1));
hold on; plot(freqs,mean(mlT_M1));
xlim([0 15]); %set(gca,'YScale','log');

subplot(1,2,2); plot(freqs,mean(meT_Cb));
hold on; plot(freqs,mean(mlT_Cb));
xlim([0 15]); %set(gca,'YScale','log');

%% Bar plot of power
figure('Color','white','Position',[404 276 1439 613]);
freq_band = {[0.1 4] [3 6] [6 14]};

for i =1:length(freq_band)
% eP_M1 = mean(mean(meT_M1(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2));
% lP_M1 = mean(mean(mlT_M1(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2));
% 
% subplot(i,2,1); bar([eP_M1,lP_M1]); hold on;
% errorbar(1,eP_M1,std(mean(meT_M1(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2))/(size(meT_M1,1)-1),'.');
% errorbar(2,lP_M1,std(mean(mlT_M1(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2))/(size(mlT_M1,1)-1),'.');

% eP_Cb = mean(mean(meT_Cb(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2));
% lP_Cb = mean(mean(mlT_Cb(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2));

eP_Cb = mean(mean(meT_Cb(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2));
lP_Cb = mean(mean(mlT_Cb(8:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2));

% subplot(i,2,2); 
bar([eP_Cb,lP_Cb]); hold on;
errorbar(1,eP_Cb,std(mean(meT_Cb(1:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2))/(size(meT_Cb,1)-1),'.');
errorbar(2,lP_Cb,std(mean(mlT_Cb(1:end,freqs>=freq_band{i}(1)&freqs<=freq_band{i}(2)),2))/(size(mlT_Cb,1)-1),'.');

end

%% Stats
x = mean(meT_Cb(:,freqs>=freq_band(1)&freqs<=freq_band(2)),2);
y = mean(mlT_Cb(:,freqs>=freq_band(1)&freqs<=freq_band(2)),2);




%% Plot all sessions

path = 'Z:\Rohit\BMI_Data\';
% sub = {'I096','I107','I110','I111','I112'};
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};
% Fs = 1.017252624511719e+03;   
% disp('running...');
% rootpath = 'Z:\Rohit\BMI_Data\';
% 
% first = [4,3,4,4,3];
% last = [6,6,10,10,10];
% robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]};
meT_all = [];
mlT_all = [];
mean_perf_consolidated = [];
for s=stroke%length(sub)
    s
%     bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
      for i=1:4
      if(s<11)
        load([path, sub{s},'\LFP_power_consolidated',num2str(i),'.mat']);
      else
         load([path, sub{s},'\Data\','\LFP_power_consolidated',num2str(i),'.mat']); 
      end
      meT_all{s,i} = meT_Cb;
      mlT_all{s,i} = mlT_Cb;
      end
end
%%
clear comp
for i =1:4
mean_perf_consolidated = [meT_all(:,i),mlT_all(:,i)];
mean_perf_consolidated = cell2mat(mean_perf_consolidated);

comp{i} = mean_perf_consolidated(:,2) - mean_perf_consolidated(:,1);
[h(i),p_val(i)] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','left');
end
%%
comp = cell2mat(comp);
[p,t,stats] = anova1(comp);
[c,m,h,gnames] = multcompare(stats); 
%%
figure('Color','white'); hold all;
bar(nanmean(mean_perf_consolidated)); 
errorbar([1,2],nanmean(mean_perf_consolidated),nanstd(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)));
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(mean_perf_consolidated','k');
%%

%% Plot all sessions - M1

path = 'Z:\Rohit\BMI_Data\';
% sub = {'I096','I107','I110','I111','I112'};
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};
% Fs = 1.017252624511719e+03;   
% disp('running...');
% rootpath = 'Z:\Rohit\BMI_Data\';
% 
% first = [4,3,4,4,3];
% last = [6,6,10,10,10];
% robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]};

M1 = [1,2,6,7,8];
meT_all = [];
mlT_all = [];
mean_perf_consolidated = [];
for s=M1%length(sub)
    s
%     bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));

      for i=1:4
      if(s<11)
        load([path, sub{s},'\LFP_power_consolidated',num2str(i),'.mat']);
      else
         load([path, sub{s},'\Data\','\LFP_power_consolidated',num2str(i),'.mat']); 
      end
      meT_all{s,i} = meT_M1;
      mlT_all{s,i} = mlT_M1;
      end

end    
%%
clear comp
for i =1:4
mean_perf_consolidated = [meT_all(:,i),mlT_all(:,i)];
mean_perf_consolidated = cell2mat(mean_perf_consolidated);
[h(i),p_val(i)] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','left');

comp{:,i} = mean_perf_consolidated(:,2) - mean_perf_consolidated(:,1); 
% Barplot
figure('Color','white'); hold all;
bar(nanmean(mean_perf_consolidated)); 
errorbar([1,2],nanmean(mean_perf_consolidated),nanstd(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)));
plot(mean_perf_consolidated','k');
end
%%
comp =cell2mat(comp);
[p,t,stats] = anova1(comp);
[c,m,h,gnames] = multcompare(stats); 
%%
% [p,t,stats] = anova1(mean_perf_consolidated);
% [c,m,h,gnames] = multcompare(stats); 
figure('Color','white'); hold all;
bar(nanmean(mean_perf_consolidated)); 
errorbar([1,2],nanmean(mean_perf_consolidated),nanstd(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)));
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(mean_perf_consolidated','k');

%% Plot all sessions - FFC

path = 'Z:\Rohit\BMI_Data\';
% sub = {'I096','I107','I110','I111','I112'};
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};
% Fs = 1.017252624511719e+03;   
% disp('running...');
% rootpath = 'Z:\Rohit\BMI_Data\';
% 
% first = [4,3,4,4,3];
% last = [6,6,10,10,10];
% robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]};

M1 = [1,2,6,7,8];
meT_all = [];
mlT_all = [];
mean_perf_consolidated = [];
for s=M1%M1%length(sub)
    s
%     bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
      clear meT mlT;
      for i=1:4
      if(s<11)
        load([path, sub{s},'\FFC_power_consolidated',num2str(i),'.mat']);
      else
         load([path, sub{s},'\Data\','\FFC_power_consolidated',num2str(i),'.mat']); 
      end

      meT_all{s,i} = meT;
      mlT_all{s,i} = mlT;
      end
end    
%%
for i =1:4
mean_perf_consolidated = [meT_all(:,i),mlT_all(:,i)];
mean_perf_consolidated = cell2mat(mean_perf_consolidated);
[h(i),p_val(i)] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','left');

comp{:,i} = mean_perf_consolidated(:,2) - mean_perf_consolidated(:,1); 
% Barplot
figure('Color','white'); hold all;
bar(nanmean(mean_perf_consolidated)); 
errorbar([1,2],nanmean(mean_perf_consolidated),nanstd(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)));
plot(mean_perf_consolidated','k');

end
%%
comp =cell2mat(comp);
[p,t,stats] = anova1(comp);
[c,m,h,gnames] = multcompare(stats); 
%%
figure('Color','white'); hold all;
bar(nanmean(mean_perf_consolidated)); 
errorbar([1,2],nanmean(mean_perf_consolidated),nanstd(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)));
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(mean_perf_consolidated','k');
%% Compare difference between early to late for intact vs stroke models

%% Plot all sessions

path = 'Z:\Rohit\BMI_Data\';
% sub = {'I096','I107','I110','I111','I112'};
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};

stroke= [4:5,9:10];

intact = [1:10];
intact = setdiff(intact, stroke);

% Fs = 1.017252624511719e+03;   
% disp('running...');
% rootpath = 'Z:\Rohit\BMI_Data\';
% 
% first = [4,3,4,4,3];
% last = [6,6,10,10,10];
% robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]};
healthy_all = [];
stroke_all = [];
mean_perf_consolidated = [];

for s=intact%length(sub)
    s
%     bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));

      if(s<11)
        load([path, sub{s},'\LFP_power_consolidated2.mat']);
      else
         load([path, sub{s},'\Data\','\LFP_power_consolidated2.mat']); 
      end
      healthy_all{s} = (mlT_Cb- meT_Cb);
%       healthy_all{s} = mlT_Cb;
%       mlT_all{s} = mlT_Cb;

end 

for s=stroke%length(sub)
    s
%     bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));


      if(s<11)
        load([path, sub{s},'\LFP_power_consolidated2.mat']);
      else
         load([path, sub{s},'\Data\','\LFP_power_consolidated2.mat']); 
      end
      stroke_all{s} = (mlT_Cb- meT_Cb);
%       stroke_all{s} = mlT_Cb;
%       mlT_all{s} = mlT_Cb;

end    

val1 = cell2mat(healthy_all(:));
val1(isnan(val1)) = [];
val1(isinf(val1)) = [];
val2 = cell2mat(stroke_all(:));
val2(isnan(val2)) = [];
val2(isinf(val2)) = [];

p_val = ranksum(val1,val2,'tail','right');
[h,p] = ttest2(val1,val2,'tail','right');
%%
% figure('Color','white'); hold all;
bar([nanmean(val1),nanmean(val2)]); hold on;
errorbar([1,2],[nanmean(val1),nanmean(val2)],[nanstd(val1)/sqrt((size(val2,1)-1)),nanstd(val2)/sqrt((size(val2,1)-1))]);
% mean_perf_consolidated = [healthy_all;stroke_all]';
% mean_perf_consolidated = cell2mat(mean_perf_consolidated);
% [h,p_val] = ttest2(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','left');
[h,p] = ttest2(val1,val2,'tail','right')

%% Anova

group = [repelem(1,length(val1)),repelem(2,length(val2))];
[p,t,stats] = anovan([(val1)',(val2)'],{group});
% [c,m,h,gnames] = multcompare(stats);
%% Plot all sessions

path = 'Z:\Rohit\BMI_Data\';
sub = {'I096','I107','I110','I111','I112'};

% Fs = 1.017252624511719e+03;   
% disp('running...');
% rootpath = 'Z:\Rohit\BMI_Data\';
% 
% first = [4,3,4,4,3];
% last = [6,6,10,10,10];
% robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]};
healthy_all = [];
M1_all = [];
mean_perf_consolidated = [];

for s=1:2%length(sub)
    s
%     bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));


      load([path, sub{s},'\LFP_power_consolidated2.mat']);
      healthy_all{s} = mlT_Cb- meT_Cb;
      M1_all{s} = mlT_M1-meT_M1;

end  

val1 = cell2mat(healthy_all(:));
val2 = cell2mat(M1_all(:));
p_val = ranksum(val1,val2,'tail','right');
%%
figure('Color','white'); hold all;
bar([mean(val1),mean(val2)]);
errorbar([1,2],[mean(val1),mean(val2)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1))],'.');
% mean_perf_consolidated = [healthy_all;stroke_all]';
% mean_perf_consolidated = cell2mat(mean_perf_consolidated);
% [h,p_val] = ttest2(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','left');

