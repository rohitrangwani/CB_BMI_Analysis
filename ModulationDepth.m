%% Get modulation depth of PSTHs from early to late trials
%  Author - Aamir Abbasi
%  Modified:  Rohit, for Cb vs M1 modulation depth comparison
%  -----------------------------------------------------------------------
%% Modulation depth of PSTHs of M1 direct units
clear;clc;close all;
disp('running');

MDdirect = [];
% disp('running...');


sessions = ["I060","I086","I096","I107"];

for s=1:length(sessions)
rootpath = char(append('Z:\Rohit\BMI_Data\M1_Data\',sessions(s),'\'));
% rootpath = 'Z:\Rohit\BMI_Data\I096\Data\';
disp(rootpath);
bmiBlocks = dir(rootpath);    

  for i=3:length(bmiBlocks) %
    
  disp(bmiBlocks(i).name);
    load([rootpath,bmiBlocks(i).name,'\Direct_Units.mat'],'PSTHdirect_tp','PSTHdirect_tn');
    psth = PSTHdirect_tp(~cellfun('isempty',PSTHdirect_tp));
    psth = [psth;PSTHdirect_tn(~cellfun('isempty',PSTHdirect_tn))];    
    for j=1:size(psth,1)
        bfr = zscore(psth{j});
        
        e_psth = bfr(1:round(size(bfr,1)/3),:);
        e_mean_psth = mean(e_psth);
        baseline_mean_psth_early = mean(e_mean_psth(1:1000));
        baseline_std_psth_early  = std(e_mean_psth(1:1000));
        peak_mean_psth_early_mx  = max(e_mean_psth(2001:6000));
        peak_mean_psth_early_mn  = min(e_mean_psth(2001:6000));
        if peak_mean_psth_early_mx > abs(peak_mean_psth_early_mn)
            e_perct_mod = round((peak_mean_psth_early_mx-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        else
            e_perct_mod = round((peak_mean_psth_early_mn-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        end
        
        l_psth = bfr(end-round(size(bfr,1)/3)+1:end,:);
        l_mean_psth = mean(l_psth);
        baseline_mean_psth_late = mean(l_mean_psth(1:1000));
        baseline_std_psth_late  = std(l_mean_psth(1:1000));
        peak_mean_psth_late_mx  = max(l_mean_psth(2001:6000));
        peak_mean_psth_late_mn  = min(l_mean_psth(2001:6000));  
        if peak_mean_psth_late_mx > abs(peak_mean_psth_late_mn)
            l_perct_mod = round((peak_mean_psth_late_mx-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        else
            l_perct_mod = round((peak_mean_psth_late_mn-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        end
        
        if peak_mean_psth_late_mx > abs(baseline_mean_psth_late+2.5*baseline_std_psth_late)...
                && peak_mean_psth_late_mx > abs(peak_mean_psth_late_mn)
            MDdirect = [MDdirect; abs(e_perct_mod) abs(l_perct_mod) 2 i]; % direct units
        end
        if peak_mean_psth_late_mn > abs(baseline_mean_psth_late-2.5*baseline_std_psth_late) ...
                && peak_mean_psth_late_mx < abs(peak_mean_psth_late_mn)
            MDdirect = [MDdirect; abs(e_perct_mod) abs(l_perct_mod) 2 i]; % direct units
        end
    end
end

a = MDdirect;
b = a(any(a,2),:);
c = b(sum(isnan(b),2)==0,:);
MDdirect = c(sum(isinf(c),2)==0,:);

% XX = [];
% for g=1:max(MDdirect(:,4))
%     a = MDdirect(MDdirect(:,4)==g,:);
%     b = a(a(:,2)>a(:,1),:);
%     if ~isempty(b)
%         XX = [XX; b(:,1) b(:,2)];
%     else
%         XX = [XX; a(:,1) a(:,2)];        
%     end
% end
% 
% [h1,b1] = hist(XX(:,1),10);
% h1 = h1/sum(h1);
% [f1,x1] = ecdf(h1);
% [h2,b2] = hist(XX(:,2),10);
% h2 = h2/sum(h2);
% [f2,x2] = ecdf(h2);

XX = [];
for g=1:max(MDdirect(:,4))
    a = MDdirect(MDdirect(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; mean(b(:,1)) mean(b(:,2))];
    else
        XX = [XX; mean(a(:,1)) mean(a(:,2))];        
    end
end

figure('Color','w','Position',[680 425 304 553]);
bar(nanmean(XX)); hold on;
errorbar(1,mean(XX(:,1)),std(XX(:,1))/sqrt(length(XX(:,1))-1));
errorbar(2,mean(XX(:,2)),std(XX(:,2))/sqrt(length(XX(:,2))-1));
plot(XX');

% e_MDdirect = splitapply(@mean,MDdirect(:,1),MDdirect(:,4));
% l_MDdirect = splitapply(@mean,MDdirect(:,2),MDdirect(:,4));
% 
% % MDdirect = MDdirect(MDdirect(:,2)-MDdirect(:,1)>0,:); 
% 
% figure;
% bar([mean(e_MDdirect) mean(l_MDdirect)]); hold on;
% errorbar(1,mean(e_MDdirect),std(e_MDdirect)/sqrt(length(e_MDdirect)-1));
% errorbar(2,mean(l_MDdirect),std(l_MDdirect)/sqrt(length(l_MDdirect)-1));
% plot([e_MDdirect l_MDdirect]')
end
%% Modulation depth PSTHs of M1 indirect units
clc;
disp('running');

MDindirect_M1 = [];
sessions = ["I060","I086","I096","I107"];

for s=1:length(sessions)
rootpath = char(append('Z:\Rohit\BMI_Data\M1_Data\',sessions(s),'\'));
% rootpath = 'Z:\Rohit\BMI_Data\I096\Data\';
disp(rootpath);
bmiBlocks = dir(rootpath);    

  for i=3:length(bmiBlocks) %
    
  disp(bmiBlocks(i).name);
    
    disp(bmiBlocks(i).name);
    load([rootpath,bmiBlocks(i).name,'\Indirect_Units_M1.mat'],'PSTHindirect');
    psth = PSTHindirect(~cellfun('isempty',PSTHindirect));
    for j=1:size(psth,1)
        bfr = zscore(psth{j});
        
        e_psth = bfr(1:round(size(bfr,1)/3),:);
        e_mean_psth = mean(e_psth);
        baseline_mean_psth_early = mean(e_mean_psth(1:1000));
        baseline_std_psth_early  = std(e_mean_psth(1:1000));
        peak_mean_psth_early_mx  = max(e_mean_psth(2001:6000));
        peak_mean_psth_early_mn  = min(e_mean_psth(2001:6000));  
        if peak_mean_psth_early_mx > abs(peak_mean_psth_early_mn)
            e_perct_mod = round((peak_mean_psth_early_mx-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        else
            e_perct_mod = round((peak_mean_psth_early_mn-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        end
        
        l_psth = bfr(end-round(size(bfr,1)/3)+1:end,:);
        l_mean_psth = mean(l_psth);
        baseline_mean_psth_late = mean(l_mean_psth(1:1000));
        baseline_std_psth_late  = std(l_mean_psth(1:1000));
        peak_mean_psth_late_mx  = max(l_mean_psth(2001:6000));
        peak_mean_psth_late_mn  = min(l_mean_psth(2001:6000));  
        if peak_mean_psth_late_mx > abs(peak_mean_psth_late_mn)
            l_perct_mod = round((peak_mean_psth_late_mx-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        else
            l_perct_mod = round((peak_mean_psth_late_mn-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        end
        
        if peak_mean_psth_late_mx > (baseline_mean_psth_late+3*baseline_std_psth_late)...
                && peak_mean_psth_late_mx > abs(peak_mean_psth_late_mn)
            MDindirect_M1 = [MDindirect_M1;abs(e_perct_mod) abs(l_perct_mod) 1, i]; % indirect units
        elseif abs(peak_mean_psth_late_mn) > abs(baseline_mean_psth_late-3*baseline_std_psth_late)...
                && peak_mean_psth_late_mx < abs(peak_mean_psth_late_mn)
            MDindirect_M1 = [MDindirect_M1;abs(e_perct_mod) abs(l_perct_mod) 1, i]; % indirect units
        else
            MDindirect_M1 = [MDindirect_M1;e_perct_mod l_perct_mod 0, i]; % unrelated units
        end
    end
end

a = MDindirect_M1;
b = a(any(a,2),:);
c = b(sum(isnan(b),2)==0,:);
MDindirect_M1 = c(sum(isinf(c),2)==0,:);

XX = [];
for g=1:max(MDindirect_M1(:,4))
    a = MDindirect_M1(MDindirect_M1(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; b(:,1) b(:,2)];
    else
        XX = [XX; a(:,1) a(:,2)];        
    end
end

[h1,b1] = hist(XX(:,1),100);
% h1 = h1/sum(h1);
[f1,x1] = ecdf(h1);
[h2,b2] = hist(XX(:,2),100);
% h2 = h2/sum(h2);
[f2,x2] = ecdf(h2);

end

%% Pie chart for direct/indirect/unrelated M1 units
XX = [];
for g=1:max(MDdirect(:,4))
    a = MDdirect(MDdirect(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; b];
    else
        XX = [XX; a];        
    end
end
for g=1:max(MDindirect_M1(:,4))
    a = MDindirect_M1(MDindirect_M1(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; b];
    else
        XX = [XX; b];        
    end
end

figure;
pie([sum(XX(:,3)==2)/size(XX,1)*100,sum(XX(:,3)==1)/size(XX,1)*100 ...
     sum(XX(:,3)==0)/size(XX,1)*100]);

 %% Pie chart for indirect/unrelated M1 units
XX = [];
for g=1:max(MDindirect_M1(:,4))
    a = MDindirect_M1(MDindirect_M1(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; b];
    else
        XX = [XX; b];        
    end
end

figure;
pie([sum(XX(:,3)==2)/size(XX,1)*100,sum(XX(:,3)==1)/size(XX,1)*100 ...
     sum(XX(:,3)==0)/size(XX,1)*100]);
%% Modulation depth PSTHs of Cb indirect units
clear;clc;close all;
disp('running');

MDindirect_Cb = [];
sessions = ["I096","I107","I110"];

blocks = {[1,2];[1,2,3];[2,6,7,10,11]};

for s=1:length(sessions)

rootpath = char(append('Z:\Rohit\BMI_Data\',sessions(s),'\Data\'));
disp(rootpath);
bmiBlocks = dir(rootpath);    

    failed_tr = zeros(1,length(bmiBlocks)-2);
    late_reward_time = zeros(1,length(bmiBlocks)-2);
    early_reward_time = zeros(1,length(bmiBlocks)-2);
  for i=1:length(blocks{s}) %
        
      b_num = 3+blocks{s}(i);    
      % Display current block id 
      disp(bmiBlocks(b_num).name);
    load([rootpath,bmiBlocks(b_num).name,'\Indirect_Units_Cb.mat'],'PSTHindirect');
    psth = PSTHindirect(~cellfun('isempty',PSTHindirect));
    for j=1:size(psth,1)
        bfr = psth{j};
        
        e_psth = bfr(1:round(size(bfr,1)/3),:);
        e_mean_psth = mean(e_psth);
        baseline_mean_psth_early = mean(e_mean_psth(1:1000));
        baseline_std_psth_early  = std(e_mean_psth(1:1000));
        peak_mean_psth_early_mx  = max(e_mean_psth(2001:4000));
        peak_mean_psth_early_mn  = min(e_mean_psth(2001:4000));  
        if peak_mean_psth_early_mx > abs(peak_mean_psth_early_mn)
            e_perct_mod = round((peak_mean_psth_early_mx-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        else
            e_perct_mod = round((peak_mean_psth_early_mn-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        end
        
        l_psth = bfr(end-round(size(bfr,1)/3)+1:end,:);
        l_mean_psth = mean(l_psth);
        baseline_mean_psth_late = mean(l_mean_psth(1:1000));
        baseline_std_psth_late  = std(l_mean_psth(1:1000));
        peak_mean_psth_late_mx  = max(l_mean_psth(2001:4000));
        peak_mean_psth_late_mn  = min(l_mean_psth(2001:4000));  
        if peak_mean_psth_late_mx > abs(peak_mean_psth_late_mn)
            l_perct_mod = round((peak_mean_psth_late_mx-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        else
            l_perct_mod = round((peak_mean_psth_late_mn-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        end
        
        if peak_mean_psth_late_mx > (baseline_mean_psth_late+2.5*baseline_std_psth_late)
            MDindirect_Cb = [MDindirect_Cb; e_perct_mod l_perct_mod 1, i]; % indirect units
        elseif peak_mean_psth_late_mn < (baseline_mean_psth_late-2.5*baseline_std_psth_late)
            MDindirect_Cb = [MDindirect_Cb; e_perct_mod l_perct_mod 1, i]; % indirect units
        else
            MDindirect_Cb = [MDindirect_Cb; e_perct_mod l_perct_mod 0, i]; % unrelated units
        end
    end
end

a = MDindirect_Cb;
b = a(any(a,2),:);
c = b(sum(isnan(b),2)==0,:);
MDindirect_Cb = c(sum(isinf(c),2)==0,:);

XX = [];
for g=1:max(MDindirect_Cb(:,4))
    a = MDindirect_Cb(MDindirect_Cb(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; mean(b(:,1)) mean(b(:,2))];
    else
        XX = [XX; mean(a(:,1)) mean(a(:,2))];        
    end
end

figure('Color','w','Position',[680 425 304 553]);
bar(mean(XX)); hold on;
errorbar(1,mean(XX(:,1)),std(XX(:,1))/sqrt(length(XX(:,1))-1));
errorbar(2,mean(XX(:,2)),std(XX(:,2))/sqrt(length(XX(:,2))-1));
plot(XX');
end



%% Modulation depth PSTHs of Cb direct units

disp('running');

MDdirect_Cb = [];
sessions = ["I096","I107","I110"];

blocks = {[1,2];[1,2,3];[2,6,7,10,11]};

for s=1:length(sessions)

rootpath = char(append('Z:\Rohit\BMI_Data\',sessions(s),'\Data\'));
disp(rootpath);
bmiBlocks = dir(rootpath);    

    failed_tr = zeros(1,length(bmiBlocks)-2);
    late_reward_time = zeros(1,length(bmiBlocks)-2);
    early_reward_time = zeros(1,length(bmiBlocks)-2);
  for i=1:length(blocks{s}) %
        
      b_num = 3+blocks{s}(i);    
      % Display current block id 
      disp(bmiBlocks(b_num).name);
    load([rootpath,bmiBlocks(b_num).name,'\Direct_units.mat'],'PSTHdirect_tp','PSTHdirect_tn');
    psth = PSTHdirect_tp(~cellfun('isempty',PSTHdirect_tp));
    psth = [psth;PSTHdirect_tn(~cellfun('isempty',PSTHdirect_tn))];    
    for j=1:size(psth,1)
        bfr = psth{j};
        
        e_psth = bfr(1:round(size(bfr,1)/3),:);
        e_mean_psth = mean(e_psth);
        baseline_mean_psth_early = mean(e_mean_psth(1:1000));
        baseline_std_psth_early  = std(e_mean_psth(1:1000));
        peak_mean_psth_early_mx  = max(e_mean_psth(2001:4000));
        peak_mean_psth_early_mn  = min(e_mean_psth(2001:4000));  
        if peak_mean_psth_early_mx > abs(peak_mean_psth_early_mn)
            e_perct_mod = round((peak_mean_psth_early_mx-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        else
            e_perct_mod = round((peak_mean_psth_early_mn-baseline_mean_psth_early)...
                /baseline_mean_psth_early*100);
        end
        
        l_psth = bfr(end-round(size(bfr,1)/3)+1:end,:);
        l_mean_psth = mean(l_psth);
        baseline_mean_psth_late = mean(l_mean_psth(1:1000));
        baseline_std_psth_late  = std(l_mean_psth(1:1000));
        peak_mean_psth_late_mx  = max(l_mean_psth(2001:4000));
        peak_mean_psth_late_mn  = min(l_mean_psth(2001:4000));  
        if peak_mean_psth_late_mx > abs(peak_mean_psth_late_mn)
            l_perct_mod = round((peak_mean_psth_late_mx-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        else
            l_perct_mod = round((peak_mean_psth_late_mn-baseline_mean_psth_late)...
                /baseline_mean_psth_late*100);
        end
        
        if peak_mean_psth_late_mx > (baseline_mean_psth_late+2.5*baseline_std_psth_late)
            MDdirect_Cb = [MDdirect_Cb; e_perct_mod l_perct_mod 1, i]; % indirect units
        elseif peak_mean_psth_late_mn < (baseline_mean_psth_late-2.5*baseline_std_psth_late)
            MDdirect_Cb = [MDdirect_Cb; e_perct_mod l_perct_mod 1, i]; % indirect units
        else
            MDdirect_Cb = [MDdirect_Cb; e_perct_mod l_perct_mod 0, i]; % unrelated units
        end
    end
end

a = MDdirect_Cb;
b = a(any(a,2),:);
c = b(sum(isnan(b),2)==0,:);
MDdirect_Cb = c(sum(isinf(c),2)==0,:);

XX = [];
for g=1:max(MDdirect_Cb(:,4))
    a = MDdirect_Cb(MDdirect_Cb(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; mean(b(:,1)) mean(b(:,2))];
    else
        XX = [XX; mean(a(:,1)) mean(a(:,2))];        
    end
end

figure('Color','w','Position',[680 425 304 553]);
bar(mean(XX)); hold on;
errorbar(1,mean(XX(:,1)),std(XX(:,1))/sqrt(length(XX(:,1))-1));
errorbar(2,mean(XX(:,2)),std(XX(:,2))/sqrt(length(XX(:,2))-1));
plot(XX');
end

%% Pie chart for direct/indirect/unrelated Cb units
XX = [];
for g=1:max(MDindirect_Cb(:,4))
    a = MDindirect_Cb(MDindirect_Cb(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; b];
    else
        XX = [XX; b];        
    end
end

figure;
pie([sum(XX(:,3)==2)/size(XX,1)*100,sum(XX(:,3)==1)/size(XX,1)*100 ...
     sum(XX(:,3)==0)/size(XX,1)*100]);
 
 %% Pie chart for direct/indirect/unrelated Cb units
XX = [];
for g=1:max(MDdirect_Cb(:,4))
    a = MDdirect_Cb(MDdirect_Cb(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; b];
    else
        XX = [XX; a];        
    end
end
for g=1:max(MDindirect_Cb(:,4))
    a = MDindirect_Cb(MDindirect_Cb(:,4)==g,:);
    b = a(a(:,2)-a(:,1)>0,:);
    if ~isempty(b)
        XX = [XX; b];
    else
        XX = [XX; b];        
    end
end

figure;
pie([sum(XX(:,3)==2)/size(XX,1)*100,sum(XX(:,3)==1)/size(XX,1)*100 ...
     sum(XX(:,3)==0)/size(XX,1)*100]);