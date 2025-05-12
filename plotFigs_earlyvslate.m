%% Rohit -2023
%% Generate performance early late bar plots for robust learning sessions - Intact

clc;clear;close all;

%Load BMI session info
bmi_session_info;

mean_perf_consolidated = [];
for s=intact%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      mean_perf_consolidated{s,n} = mean_perf(1:2);
    end
end
perf_1 = mean_perf_consolidated;
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated1 = mean_perf_consolidated;
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right');
%% plot and compare
[p,t,stats] = anova1(mean_perf_consolidated);
[c,m,h,gnames] = multcompare(stats); 
figure('Color','white'); hold all;
bar(mean(mean_perf_consolidated),0.6,'w');
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
%% linear mixed effect model analysis - intact

rat_id = [];
for col=1:(size(intersection,2))
    if( ~isempty(find(col==intact)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end

%
tbl = [[mean_perf_consolidated1(:,1);mean_perf_consolidated1(:,2)],[zeros(length(mean_perf_consolidated1),1);...
    ones(length(mean_perf_consolidated1),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)


%% Generate performance early late bar plots for robust learning sessions - Stroke

close all;

%Load BMI session info
bmi_session_info;

mean_perf_consolidated = [];
for s=stroke
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      mean_perf_consolidated{s,n} = mean_perf(1:2);
    end
end
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right');

mean_perf_consolidated2 = mean_perf_consolidated;
%%
[p,t,stats] = anova1(mean_perf_consolidated);
[c,m,h,gnames] = multcompare(stats); 
figure('Color','white'); hold all;
bar(mean(mean_perf_consolidated),0.6,'w');
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
%% linear mixed effect model analysis - stroke

rat_id = [];
for col=1:(size(intersection,2))
    if( ~isempty(find(col==stroke)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end

%
tbl = [[mean_perf_consolidated2(:,1);mean_perf_consolidated2(:,2)],[zeros(length(mean_perf_consolidated2),1);...
    ones(length(mean_perf_consolidated2),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)

%% Generate performance early late bar plots for robust learning sessions - M1

close all;
 
% sub = {'I060','I061','I086','I096','I1107'};

Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\M1_Data\';

mean_perf_consolidated = [];
% for s=1:length(sub)
%     s
%     bmiBlocks = dir(strcat(rootpath,sub{s},'\'));
%     for n=3:length(bmiBlocks)
% 
%       load([rootpath, sub{s},'\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
%       mean_perf_consolidated{s,n} = mean_perf(1:2);
%     end
% end
T = readmatrix([rootpath,'BMI_perf.xlsx']);
% mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated = T(:,1:2);
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right');
mean_perf_consolidated3 = mean_perf_consolidated;
%%
[p,t,stats] = anova1(mean_perf_consolidated);
[c,m,h,gnames] = multcompare(stats); 
figure('Color','white'); hold all;
bar(mean(mean_perf_consolidated),0.6,'w');
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');

%% linear mixed effect model analysis - M1

rat_id = T(:,3);
rat_id = (rat_id - 60);
tbl = [[mean_perf_consolidated3(:,1);mean_perf_consolidated3(:,2)],[zeros(length(mean_perf_consolidated3),1);...
    ones(length(mean_perf_consolidated3),1)],[rat_id;rat_id]];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)

%% Compare all 3 , cb-healthy vs stroke vs M1

% perc1 = (mean_perf_consolidated1(:,1)-mean_perf_consolidated1(:,2))/mean_perf_consolidated1(:,1);
% perc2 = (mean_perf_consolidated2(:,1)-mean_perf_consolidated2(:,2))/mean_perf_consolidated2(:,1);
% perc3 = (mean_perf_consolidated3(:,1)-mean_perf_consolidated3(:,2))/mean_perf_consolidated3(:,1);


perc1 = mean_perf_consolidated1(:,2);
perc2 = mean_perf_consolidated2(:,2);
perc3 = mean_perf_consolidated3(:,2);

%% box plot
figure;
v = [1 2 3];
groups = repelem(v,[length(perc1)  length(perc2) length(perc3)]);
h = boxplot([perc1;perc2;perc3],groups,'symbol','');
% h = boxplot([perc1(:,1);perc2(:,24);perc3(:,12)],groups,'symbol','');
% h = boxplot([perc1(:,13);perc2(:,2);perc3(:,6)],groups,'symbol','');
set(gca,'LineWidth',2)
set(findobj(gca,'type','line'),'linew',2)
% % set(gca, 'FontSize', 14); 
% set(h(5,1), 'YData', [q10 q90 q90 q10 q10]);% blue box  
% 
% upWhisker = get(h(1,1), 'YData');  
% set(h(1,1), 'YData', [q90 upWhisker(2)])  
% 
% dwWhisker = get(h(2,1), 'YData');  
% set(h(2,1), 'YData', [ dwWhisker(1) q10])  
%%
close all;
v = [1 2 3];
groups = repelem(v,[length(perc1)  length(perc2) length(perc3)]);
% [p,t,stats] = anova1([perc1(:,1);perc2(:,24);perc3(:,12)],groups);
% [p,t,stats] = anova1([perc1(:,1);perc2(:,20);perc3(:,12)],groups);
[p,t,stats] = kruskalwallis([perc1;perc2;perc3],groups)
set(gca,'LineWidth',2)
set(gca,'TickDir','out');
set(findobj(gca,'type','line'),'linew',3)
% boxplot('symbol', '')
figure;
[c,m,h,gnames] = multcompare(stats);




%% Compare session with Tp only vs Tp and Tn

clc;clear;close all;

%Load BMI session info
bmi_session_info;

% Tp only sessions
mean_perf_consolidated = [];
for s=stroke%4:length(sub)
    
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat'],'Waves_tn');
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');

      if (sum(sum(~cellfun(@isempty,Waves_tn)))) == 0
        mean_perf_consolidated{s,n} = mean_perf(1:2);
      end
    end
end
perf_1 = mean_perf_consolidated;
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated1 = mean_perf_consolidated;
figure; bar(mean(mean_perf_consolidated),0.6,'w'); hold on;
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
hold on; plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right')

% Tp + Tn sessions
mean_perf_consolidated = [];
for s=stroke%4:length(sub)
    
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat'],'Waves_tn');
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');

     if (sum(sum(~cellfun(@isempty,Waves_tn)))) ~= 0
        mean_perf_consolidated{s,n} = mean_perf(1:2);
     end
    end
end
perf_2 = mean_perf_consolidated;
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated2 = mean_perf_consolidated;
figure; bar(mean(mean_perf_consolidated),0.6,'w'); hold on;
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
hold on; plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right')

perc1 = (mean_perf_consolidated1(:,1)-mean_perf_consolidated1(:,2))./mean_perf_consolidated1(:,1);
perc2 = (mean_perf_consolidated2(:,1)-mean_perf_consolidated2(:,2))./mean_perf_consolidated2(:,1);
v = [1 2];
groups = repelem(v,[length(perc1)  length(perc2)]);

[p,t,stats] = kruskalwallis([perc1;perc2],groups)
ylim([0,1]);
set(gca,'TickDir','out');
%%
perc1 = mean_perf_consolidated1(:,2);
perc2 = mean_perf_consolidated2(:,2);
%
close all;
v = [1 2];
groups = repelem(v,[length(perc1)  length(perc2)]);
% [p,t,stats] = anova1([perc1(:,1);perc2(:,24);perc3(:,12)],groups);
% [p,t,stats] = anova1([perc1(:,1);perc2(:,20);perc3(:,12)],groups);
[p,t,stats] = anova1([perc1;perc2],groups);
set(gca,'LineWidth',2)
set(gca,'TickDir','out');
set(findobj(gca,'type','line'),'linew',3)
% boxplot('symbol', '')
figure;
[c,m,h,gnames] = multcompare(stats);


%% Compare session with only 1 unit in Tp vs everthing else

clc;clear;close all;

%Load BMI session info
bmi_session_info;

% Tp only sessions
mean_perf_consolidated = [];
for s=stroke%4:length(sub)
    
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat'],'Waves_tp');
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');

      if (sum(sum(~cellfun(@isempty,Waves_tp)))) == 1
        mean_perf_consolidated{s,n} = mean_perf(1:2);
      end
    end
end
perf_1 = mean_perf_consolidated;
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated1 = mean_perf_consolidated;

figure; bar(mean(mean_perf_consolidated),0.6,'w');hold on;
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
hold on; plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right')

% Tp + Tn sessions
mean_perf_consolidated = [];
for s=stroke%4:length(sub)
    
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat'],'Waves_tp');
      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');

     if (sum(sum(~cellfun(@isempty,Waves_tp)))) > 1
        mean_perf_consolidated{s,n} = mean_perf(1:2);
     end
    end
end
perf_2 = mean_perf_consolidated;
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated2 = mean_perf_consolidated;

figure; bar(mean(mean_perf_consolidated),0.6,'w');hold on;
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
hold on; plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right')

%%
% perc1 = mean_perf_consolidated1(:,2);
% perc2 = mean_perf_consolidated2(:,2);

perc1 = (mean_perf_consolidated1(:,1)-mean_perf_consolidated1(:,2))./mean_perf_consolidated1(:,1);
perc2 = (mean_perf_consolidated2(:,1)-mean_perf_consolidated2(:,2))./mean_perf_consolidated2(:,1);
v = [1 2];
groups = repelem(v,[length(perc1)  length(perc2)]);

[p,t,stats] = kruskalwallis([perc1;perc2],groups)
ylim([0,1]);
set(gca,'TickDir','out');
%%
close all;
v = [1 2];
groups = repelem(v,[length(perc1)  length(perc2)]);
% [p,t,stats] = anova1([perc1(:,1);perc2(:,24);perc3(:,12)],groups);
% [p,t,stats] = anova1([perc1(:,1);perc2(:,20);perc3(:,12)],groups);
[p,t,stats] = anova1([perc1;perc2],groups);
set(gca,'LineWidth',2)
set(gca,'TickDir','out');
set(findobj(gca,'type','line'),'linew',3)
% boxplot('symbol', '')
figure;
[c,m,h,gnames] = multcompare(stats);







%% Compare purkinje direct unit sessions to all 
clc;clear;close all;

%Load BMI session info
bmi_session_info;

mean_perf_consolidated = [];
for s=[4,12]%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    if s == 4
        session = 3;
    else
        session = 9;
    end
    for n=session

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      mean_perf_consolidated{s,n} = mean_perf(1:2);
    end
end
perf_1 = mean_perf_consolidated;
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated1 = mean_perf_consolidated;
figure; bar(mean(mean_perf_consolidated),0.6,'w');hold on;
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
hold on; plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right')

mean_perf_consolidated = [];
for s=stroke%4:length(sub)
    s
    if s == 4
        session = [4,5,8,9,10];
    elseif s == 12
        session = [7,11,14];
    else
        session = intersection{s};
    end

    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));

    for n=session

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat'],'mean_perf');
      mean_perf_consolidated{s,n} = mean_perf(1:2);
    end
end
perf_2 = mean_perf_consolidated;
mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated2 = mean_perf_consolidated;
figure; bar(mean(mean_perf_consolidated),0.6,'w');hold on;
errorbar([1,2],mean(mean_perf_consolidated),std(mean_perf_consolidated)/sqrt((size(mean_perf_consolidated,1)-1)),'color','black','linewidth',4);
hold on; plot(mean_perf_consolidated','k');
ylim([0,15]);
set(gca,'TickDir','out');
[h,p_val] = ttest(mean_perf_consolidated(:,1),mean_perf_consolidated(:,2),'Tail','right')

% perc1 = mean_perf_consolidated1(:,2);
% perc2 = mean_perf_consolidated2(:,2);
perc1 = (mean_perf_consolidated1(:,1)-mean_perf_consolidated1(:,2))./mean_perf_consolidated1(:,1);
perc2 = (mean_perf_consolidated2(:,1)-mean_perf_consolidated2(:,2))./mean_perf_consolidated2(:,1);
v = [1 2];
groups = repelem(v,[length(perc1)  length(perc2)]);
%%
[p,t,stats] = kruskalwallis([perc1;perc2],groups)
ylim([0,1]);
set(gca,'TickDir','out');
%
% [h,p_val] = ttest2(mean_perf_consolidated1(:,2),mean_perf_consolidated2(:,2))
