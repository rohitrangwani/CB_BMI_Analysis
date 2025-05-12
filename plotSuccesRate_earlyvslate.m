%% plot success rate in early vs late trials for all and analyze
% Rohit

%% Get success rate - early late bar plots for robust learning sessions - Intact

clc;clear;close all;

%Load BMI session info
bmi_session_info;

tOutcome = [];
for s=intact%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat']);
      early = valid_perf(1:floor(length(valid_perf)/3));
      late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 

      tOutcome{s,n} = [sum(early<15)/sum(length(early))*100,sum(late<15)/sum(length(late))*100]; % successful early and late
%       tOutcome(s,n,2) = sum(late==15)/sum(length(late))*100; % successful late
      clear valid_pref 
    end
end
tOutcome = cell2mat(tOutcome(:));
tOutcome1 = tOutcome;
%%
figure;
[p,t,stats] = anova1(tOutcome);
[c,m,h,gnames] = multcompare(stats); 
figure('Color','white'); hold all;
bar(mean(tOutcome),0.6,'w');
errorbar([1,2],mean(tOutcome),std(tOutcome)/sqrt((size(tOutcome,1)-1)),'color','black','linewidth',4);
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(tOutcome','k');
ylim([0,100]);
xlim([0.5,2.5]);
set(gca,'TickDir','out')

%% linear mixed effect model analysis - intact

rat_id = [];
for col=1:(size(intersection,2))
    if( ~isempty(find(col==intact)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end

%
tbl = [[tOutcome1(:,1);tOutcome1(:,2)],[zeros(length(tOutcome1),1);...
    ones(length(tOutcome1),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)
%%
%% Get success rate - early late bar plots for robust learning sessions - Stroke

%Load BMI session info
bmi_session_info;

tOutcome = [];
for s=stroke%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%robust_session{s}%length(bmiBlocks)-1

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat']);
      early = valid_perf(1:floor(length(valid_perf)/3));
      late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 

      tOutcome{s,n} = [sum(early<15)/sum(length(early))*100,sum(late<15)/sum(length(late))*100]; % successful early and late
%       tOutcome(s,n,2) = sum(late==15)/sum(length(late))*100; % successful late
      clear valid_pref 
    end
end
tOutcome = cell2mat(tOutcome(:));
tOutcome2 = tOutcome;
%%
figure;

[p,t,stats] = anova1(tOutcome);
[c,m,h,gnames] = multcompare(stats); 
figure('Color','white'); hold all;
bar(mean(tOutcome),0.6,'w');
errorbar([1,2],mean(tOutcome),std(tOutcome)/sqrt((size(tOutcome,1)-1)),'color','black','linewidth',4);
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(tOutcome','k');
ylim([0,100]);
xlim([0.5,2.5]);
set(gca,'TickDir','out')
%% linear mixed effect model analysis - stroke

rat_id = [];
for col=1:(size(intersection,2))
    if( ~isempty(find(col==stroke)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end

%
tbl = [[tOutcome2(:,1);tOutcome2(:,2)],[zeros(length(tOutcome2),1);...
    ones(length(tOutcome2),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)
%% Generate performance early late bar plots for robust learning sessions - M1

close all;
 
Fs = 1.017252624511719e+03;   
disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\M1_Data\';

tOutcome = [];
T = readmatrix([rootpath,'BMI_perf.xlsx']);
tOutcome = T(:,4:5);
tOutcome =  ones(size(tOutcome))*100-tOutcome;
[h,p_val] = ttest(tOutcome(:,1),tOutcome(:,2),'Tail','right');
tOutcome3 = tOutcome;
%%
figure;

[p,t,stats] = anova1(tOutcome);
[c,m,h,gnames] = multcompare(stats); 
figure('Color','white'); hold all;
bar(mean(tOutcome),0.6,'w');
errorbar([1,2],mean(tOutcome),std(tOutcome)/sqrt((size(tOutcome,1)-1)),'color','black','linewidth',4);
% scatter(ones(size(mean_perf_consolidated,1),1),mean_perf_consolidated(:,1));
% scatter(ones(size(mean_perf_consolidated,1),1)*2,mean_perf_consolidated(:,2));
plot(tOutcome','k');
ylim([0,100]);
xlim([0.5,2.5]);
set(gca,'TickDir','out')

%% linear mixed effect model analysis - M1

rat_id = T(:,3);
rat_id = (rat_id - 60);
tbl = [[tOutcome3(:,1);tOutcome3(:,2)],[zeros(length(tOutcome3),1);...
    ones(length(tOutcome3),1)],[rat_id;rat_id]];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)


%% Compare all 3 , cb-healthy vs stroke vs M1
% 
% perc1 = (tOutcome1(:,2)-tOutcome1(:,1))/tOutcome1(:,2);
% perc2 = (tOutcome2(:,2)-tOutcome2(:,1))/tOutcome2(:,2);
% perc3 = (tOutcome3(:,2)-tOutcome3(:,1))/tOutcome3(:,2);


perc1 = tOutcome1(:,2);
perc2 = tOutcome2(:,2);
perc3 = tOutcome3(:,2);

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
% [p,t,stats] = anova1([perc1(:,1);perc2(:,1);perc3(:,1)],groups);
% [p,t,stats] = anova1([perc1(:,1);perc2(:,20);perc3(:,12)],groups);
[p,t,stats] = kruskalwallis([perc1;perc2;perc3],groups);
set(gca,'LineWidth',2)
set(gca,'TickDir','out');
set(findobj(gca,'type','line'),'linew',3)
% boxplot('symbol', '')
figure;
[c,m,h,gnames] = multcompare(stats);
