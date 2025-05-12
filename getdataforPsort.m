%%%% Author - Rohit 2023
%%%% BMI Data Analysis Gulati Lab
%%%% SCRIPT TO READ CONTINOUS SPIKING DATA FROM RS4 
%% Read TDT blocks and save channel data to a mat file for psort

% clc; clear; close;
% disp('running...');
addpath(genpath('D:\TDTDataReader'));
root = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\'; 
savepath = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\';
cd(root);

      block =  dir(root);
      blockpath = block(3).name;
      ch = 11;  
        
      % Read Cb raw data from RS4 .sev files (Cb channels are saved as ch 33-96 by RS4)
      raw_Cb = SEV2mat(blockpath,'CHANNEL', ch);
      
%       Extract Cb single units continous data
%       su_Cb = raw_Cb.RSn1.data;
      
      ch_data = raw_Cb.RSn1.data;
      
      sample_rate = raw_Cb.RSn1.fs;
      ch_time = 0:1/sample_rate:((length(ch_data)-1)/sample_rate);
      
      save(strcat(root,blockpath,'_Ch_',num2str(ch),'.mat'),'ch_time','ch_data','sample_rate');
      
      
%% Read all relevant TDT blocks and save channel data to a mat file for psort
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154'};

Cb_arr = [1,1,0,0,0,1,1,1,0,0,0]; %% Set to 1 for sessions if Cb channels are after M1 channels

robust_session = {[4],[3],[5,7,8,10,13],[3,4,5,8],[5,6,9],...
    [4],[5,7,8],[8],[7],[3,5,6],[8,14,15,16,18,19]};

p_ch = {[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],...
    [12], [8, 1], [1], [1, 14], [5, 7, 13, 21, 26], [34, 3, 7, 39, 40, 10, 12, 16, 29],...
    [1, 3, 9, 12, 16, 17, 18, 20, 23, 25, 27, 28, 29, 33, 36, 39, 41, 43, 44, 46, 53, 55, 61],...
    [35, 36, 39, 41, 43, 44, 46, 61, 21, 53, 55, 29], [32, 4, 38, 10, 42, 43, 51, 61, 21, 27, 29, 31],...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],...
    [1, 2, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24],...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20], [11, 31],...
    [3, 4, 5, 35, 37, 8, 40, 43, 12, 44, 17, 20, 54, 29], [40, 43, 44, 13, 54], ...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], [1, 2, 3, 9, 12, 13, 15],...
    [1, 33, 8, 9, 10, 11, 15, 16, 17, 18, 22, 23, 25], [33], [16, 33], ...
    [1, 2, 5, 7, 15, 19, 20, 22, 23, 25, 27], [1, 2, 3, 4, 5, 9, 10, 11, 13, 17, 18, 21], ...
    [4, 5, 6, 7, 10, 11, 14, 16, 17], [2, 3, 4, 7, 9, 10, 11, 13, 14, 15, 16], ...
    [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 17, 21], [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 15, 18, 20, 21]};
% clc; clear; close;
% disp('running...');

count = 24;%length(cell2mat(robust_session))
for s=11:length(sub)
    s
    rootpath = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Rohit\BMI_Data\';
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}
    
% addpath(genpath('D:\TDTDataReader'));
% root = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\'; 
% savepath = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\';
% cd(root);

      root = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\';
      blockpath = strcat(root,sub{s},'\RS4_Data\',bmiBlocks(n).name);
      if Cb_arr(s) == 1
           p_ch{count} = p_ch{count} +32;
           ref = SEV2mat(blockpath,'CHANNEL',(33:96));
           ref_data = ref.RSn1.data;
           ref_data = median(ref_data);
      else
           ref = SEV2mat(blockpath,'CHANNEL',(1:64));
           ref_data = ref.RSn1.data;
           ref_data = median(ref_data);
      end

      for ch=p_ch{count}

           % Read Cb raw data from RS4 .sev files (Cb channels are saved as ch 33-96 by RS4)
          raw_Cb = SEV2mat(blockpath,'CHANNEL', ch+1);
          
    %       Extract Cb single units continous data
    %       su_Cb = raw_Cb.RSn1.data;

          ch_data = raw_Cb.RSn1.data;
          ch_data = ch_data - ref_data;
          sample_rate = raw_Cb.RSn1.fs;
          ch_time = 0:1/sample_rate:((length(ch_data)-1)/sample_rate);

          save(strcat(rootpath,'\psort_data\I154\',bmiBlocks(n).name,'_Ch_',num2str(ch),'.mat'),'ch_time','ch_data','sample_rate');
      end
      count = count +1;
    end
end


%% Read all relevant TDT blocks and save channel data to a mat file for psort
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160'};

Cb_arr = [1,1,0,0,0,1,1,1,0,0,0,0]; %% Set to 1 for sessions if Cb channels are after M1 channels

robust_session = {[4],[3],[5,7,8,10,13],[3,4,5,8],[5,6,9],...
    [4],[5,7,8],[8],[7],[3,5,6],[8,14,15,16,18,19],[7,9,11,14]};

p_ch = {[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],...
    [12], [8, 1], [1], [1, 14], [5, 7, 13, 21, 26], [34, 3, 7, 39, 40, 10, 12, 16, 29],...
    [1, 3, 9, 12, 16, 17, 18, 20, 23, 25, 27, 28, 29, 33, 36, 39, 41, 43, 44, 46, 53, 55, 61],...
    [35, 36, 39, 41, 43, 44, 46, 61, 21, 53, 55, 29], [32, 4, 38, 10, 42, 43, 51, 61, 21, 27, 29, 31],...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],...
    [1, 2, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24],...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20], [11, 31],...
    [3, 4, 5, 35, 37, 8, 40, 43, 12, 44, 17, 20, 54, 29], [40, 43, 44, 13, 54], ...
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], [1, 2, 3, 9, 12, 13, 15],...
    [1, 33, 8, 9, 10, 11, 15, 16, 17, 18, 22, 23, 25], [33], [16, 33], ...
    [1, 2, 5, 7, 15, 19, 20, 22, 23, 25, 27], [1, 2, 3, 4, 5, 9, 10, 11, 13, 17, 18, 21], ...
    [4, 5, 6, 7, 10, 11, 14, 16, 17], [2, 3, 4, 7, 9, 10, 11, 13, 14, 15, 16], ...
    [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 17, 21], [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 15, 18, 20, 21],...    
    [2, 3, 4, 5, 6, 9, 11, 13, 14, 15, 19, 22, 24, 25, 26, 27, 29, 30, 31, 33, 37, 39, 40, 41, 43, 45, 46, 49,...
    51, 52, 54, 56, 57, 58, 59, 61, 62], [2, 3, 4, 5, 6, 9, 13, 14, 15, 19, 20, 22, 24, 25, 29, 30, 31, 32, 34,...
    35, 37, 38, 39, 42, 43, 45, 46, 48, 51, 52, 54, 56, 57, 58, 60, 61, 62],...
    [1, 2, 4, 5, 6, 9, 12, 13, 14, 15, 18, 19, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 37, 38, 39,...
    40, 41, 42, 43, 45, 46, 48, 50, 51, 52, 55, 57, 58, 59, 60, 61, 62, 63],...
    [3, 4, 6, 8, 9, 11, 24, 30, 32, 35, 36, 38, 39, 40, 41, 43, 44, 45, 46, 48, 49, 50, 51, 52, 59, 61, 62, 63]};
% clc; clear; close;
% disp('running...');

count = 30;
for s=12:length(sub)
    s
    rootpath = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Rohit\BMI_Data\';
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}
    
% addpath(genpath('D:\TDTDataReader'));
% root = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\'; 
% savepath = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\';
% cd(root);

      root = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\';
      blockpath = strcat(root,sub{s},'\RS4_Data\',bmiBlocks(n).name);
      if Cb_arr(s) == 1
           p_ch{count} = p_ch{count} +32;
           ref = SEV2mat(blockpath,'CHANNEL',(33:96));
           ref_data = ref.RSn1.data;
           ref_data = median(ref_data);
      else
           ref = SEV2mat(blockpath,'CHANNEL',(1:64));
           ref_data = ref.RSn1.data;
           ref_data = median(ref_data);
      end

      for ch=p_ch{count}

           % Read Cb raw data from RS4 .sev files (Cb channels are saved as ch 33-96 by RS4)
          raw_Cb = SEV2mat(blockpath,'CHANNEL', ch+1);
          
    %       Extract Cb single units continous data
    %       su_Cb = raw_Cb.RSn1.data;

          ch_data = raw_Cb.RSn1.data;
          ch_data = ch_data - ref_data;
          sample_rate = raw_Cb.RSn1.fs;
          ch_time = 0:1/sample_rate:((length(ch_data)-1)/sample_rate);

          save(strcat(rootpath,'\psort_data\I160\',bmiBlocks(n).name,'_Ch_',num2str(ch),'.mat'),'ch_time','ch_data','sample_rate');
      end
      count = count +1;
    end
end


%% Read all relevant TDT blocks and save channel data to a mat file for psort: after I160 with updated channel ids
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160','I161','I170'};

Cb_arr = [1,1,0,0,0,1,1,1,0,0,0,0,1,1]; %% Set to 1 for sessions if Cb channels are after M1 channels

robust_session = {[4,5],[3,4,6],[7,9,10,13,14],[3,4,5,8,9,10],[5,7,9],...
    [4,6,8],[5,6,7,8,9],[8],[],[3,6],[8,14,16,18],[7,9,11,14],[8,9],[9]};

%I160 and after
p_ch = {[2, 3, 4, 5, 6, 9, 11, 13, 14, 15, 19, 22, 24, 25, 26, 27, 29, 30, 31, 33, 37, 39, 40, 41, 43, 45, 46, 49, 51, 52, 54, 56, 57, 58, 59, 61, 62],...
    [2, 3, 4, 5, 6, 9, 13, 14, 15, 19, 20, 22, 24, 25, 29, 30, 31, 32, 34, 35, 37, 38, 39, 42, 43, 45, 46, 48, 51, 52, 54, 56, 57, 58, 60, 61, 62],...
    [1, 2, 4, 5, 6, 9, 12, 13, 14, 15, 18, 19, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 37, 38, 39, 40, 41, 42, 43, 45, 46, 48, 50, 51, 52, 55, 57, 58, 59, 60, 61, 62, 63],...
    [3, 4, 6, 8, 9, 11, 24, 30, 32, 35, 36, 38, 39, 40, 41, 43, 44, 45, 46, 48, 49, 50, 51, 52, 59, 61, 62, 63],...
    [2, 5, 7, 8, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 26, 27, 29, 30, 32, 36, 40, 41, 42, 45, 46, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63],...
    [1, 2, 3, 4, 6, 7, 11, 12, 13, 14, 15, 23, 26, 29, 30, 46, 48, 50, 51, 52, 53, 54, 55, 58, 60, 62, 63],...
    [0, 1, 2, 3, 7, 8, 10, 12, 13, 15, 22, 23, 24, 27, 28, 29, 31, 35, 37, 40, 41, 44, 45, 53, 56, 57, 59, 61, 62, 63]};
% clc; clear; close;
% disp('running...');

count = 1;
for s=12:length(sub)
    s
    rootpath = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Rohit\BMI_Data\';
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}
    
% addpath(genpath('D:\TDTDataReader'));
% root = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\'; 
% savepath = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\I111\RS4_Data\';
% cd(root);


      root = '\\smb.researchcampus.csmc.edu\gulatitlab-data\Cb_BMI\';
      blockpath = strcat(root,sub{s},'\RS4_Data\',bmiBlocks(n).name);
      if Cb_arr(s) == 1
           p_ch{count} = p_ch{count} +32;
           ref = SEV2mat(blockpath,'CHANNEL',(33:96));
           ref_data = ref.RSn1.data;
           ref_data = median(ref_data);
      else
           ref = SEV2mat(blockpath,'CHANNEL',(1:64));
           ref_data = ref.RSn1.data;
           ref_data = median(ref_data);
      end

      for ch=p_ch{count}

           % Read Cb raw data from RS4 .sev files (Cb channels are saved as ch 33-96 by RS4)
          raw_Cb = SEV2mat(blockpath,'CHANNEL', ch+1);
          
    %       Extract Cb single units continous data
    %       su_Cb = raw_Cb.RSn1.data;

          ch_data = raw_Cb.RSn1.data;
          ch_data = ch_data - ref_data;
          sample_rate = raw_Cb.RSn1.fs;
          ch_time = 0:1/sample_rate:((length(ch_data)-1)/sample_rate);

          save(strcat(rootpath,'\psort_data\new\',bmiBlocks(n).name,'_Ch_',num2str(ch),'.mat'),'ch_time','ch_data','sample_rate');
      end
      count = count +1;
    end
end