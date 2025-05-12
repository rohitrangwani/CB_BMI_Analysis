%% Script to plot the sucesss rate data for reach
% Rohit - 2024

rootpath = 'Z:\Rohit\BMI_Data\Results\';
data = readmatrix([rootpath,'ReachData.xlsx']);

%% Scatter plot
% c = linspace(1,10,7);
figure;
hold all;

cmap = parula(9);
c=1:10;
% set(gca,'colororder',parula(32))

scatter(data(6:54,1),data(6:54,2:10),100,cmap,'filled'); 
colormap(gca,"winter");
%     scatter(data([2:4,6:end],1),data([2:4,6:end],i),100,'filled'); hold on;

% scatter(data(2,3:8),repelem(0.5,6),'filled'); hold on;

ylim([0,2]);
xlim([-3,40]);

yline(mean(data(6,2:10)));
set(gca, 'TickDir', 'out')

%% Get BMI performance data
%Load BMI session info
bmi_session_info;

mean_perf_consolidated = [];
tOutcome = [];

good_rec = [4,9,11,13,14,15];
bad_rec = [5,10,12];

for s=stroke
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%intersection{s}%length(bmiBlocks)-1

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat']);
      mean_perf_consolidated{s,n} = mean_perf(1:2);
      early = valid_perf(1:floor(length(valid_perf)/3));
      late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 

      tOutcome{s,n} = [sum(early<15)/sum(length(early))*100,sum(late<15)/sum(length(late))*100]; % successful early and late
      clear valid_pref mean_perf 
    end
end

mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated1 = mean_perf_consolidated;

tOutcome = cell2mat(tOutcome(:));
tOutcome1 = tOutcome;

%% get BMI days

days = data(55:60,1:10);
% days(isnan(days)) = [];
% idx  =isnan(days);
% days(idx) =[];
days = days(~isnan(days));
%% get reach data
reach_success = [];
count = 1; 

%all stroke rats
for i=2:10
    reach_success{i} =  data(days(count:count+length(robust_session{stroke(i-1)})-1)+6,i);
    days(count:count+length(robust_session{stroke(i-1)})-1)
    count = count + length(robust_session{stroke(i-1)});
end

% good_rec_col = [2,4,6,8,9,10];
%good recovery rats only
% for i=2:7
%     reach_success{i} =  data(days(count:count+length(robust_session{good_rec(i-1)})-1)+6,good_rec_col(i-1));
%     days(count:count+length(robust_session{good_rec(i-1)})-1)
%     count = count + length(robust_session{good_rec(i-1)});
% end

% bad_rec_col = [3,5,7];
% bad recovery rats only
% for i=bad_rec_col
%     reach_success{i} =  data(days(count:count+length(robust_session{bad_rec(i-1)})-1)+6,i);
%     days(count:count+length(robust_session{bad_rec(i-1)})-1)
%     count = count + length(robust_session{bad_rec(i-1)});
% end

% days = reshape(days,1,size(days,1)*size(days,2));
reach_success = cell2mat(reach_success(:));
% Correlation of reach success with BMI performance

bad_reach = ~good_reach;
perc1 = (mean_perf_consolidated1(:,1)-mean_perf_consolidated1(:,2))/mean_perf_consolidated1(:,1);
perc2 = (tOutcome1(:,2)-tOutcome1(:,1))/tOutcome1(:,2);

perc3 = perc1(:,8);%perc1(:,4);%perc1(:,5);%perc1(:,8);
perc4 = perc2(:,1);%perc2(:,6);%perc2(:,1);
idx  = isnan(reach_success);
reach_success(idx) = [];

perc3(idx) = [];
perc4(idx) = [];

figure();
scatter(reach_success,perc3,'k','filled')
[r,p] = corrcoef(reach_success,perc3)
set(gca, 'TickDir', 'out')
figure();
scatter(reach_success,perc4,'k','filled')
[r,p] = corrcoef(reach_success,perc4)
set(gca, 'TickDir', 'out')


%% good reachers only
good_reach = [1;1;1;1;0;1;1;1;1;0;0;1;1;1;1;1];

g_reach_success = reach_success(logical(good_reach'));
g_perc3 = perc3(logical(good_reach'));
g_perc4 = perc4(logical(good_reach'));

figure();
scatter(g_perc3,g_reach_success,'k','filled')
[r,p] = corrcoef(g_reach_success,g_perc3)
set(gca, 'TickDir', 'out')
xlim([0,0.4])
figure();
scatter(g_perc4,g_reach_success,'k','filled')
[r,p] = corrcoef(g_reach_success,g_perc4)
set(gca, 'TickDir', 'out')
xlim([0,0.4])

%% bad reachers only

good_reach = [1;1;1;1;0;1;1;1;1;0;0;1;1;1;1;1];
bad_reach = ~good_reach;

b_reach_success = reach_success(logical(bad_reach'));
b_perc3 = perc3(logical(bad_reach'));
b_perc4 = perc4(logical(bad_reach'));

figure();
scatter(b_perc3,b_reach_success,'k','filled')
[r,p] = corrcoef(b_reach_success,b_perc3)
set(gca, 'TickDir', 'out')
xlim([0,0.4])
ylim([0,0.7])
figure();
scatter(b_perc4,b_reach_success,'k','filled')
[r,p] = corrcoef(b_reach_success,b_perc4)
set(gca, 'TickDir', 'out')
xlim([0,0.4])
ylim([0,0.7])

%% Linear mixed effect model
rat_id = [1;1;1;1;2;3;3;5;5;6;6;7;8;8;9;9];
tbl = [perc3,perc4,...
   reach_success,good_reach,rat_id];
tbl = array2table(tbl,'VariableNames',{'time','success','reach','stroke_rec','RatID'});
formula = 'reach ~ time + success + stroke_rec + (time + success + stroke_rec | RatID)';
lme_R2 = fitlme(tbl,formula)
%% Get BMI performance data
%Load BMI session info
bmi_session_info;

mean_perf_consolidated = [];
tOutcome = [];

good_rec = [4,9,11,13,14,15];
bad_rec = [5,10,12];

for s=good_rec%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat']);
      mean_perf_consolidated{s,n} = mean_perf(1:2);
      early = valid_perf(1:floor(length(valid_perf)/3));
      late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 

      tOutcome{s,n} = [sum(early<15)/sum(length(early))*100,sum(late<15)/sum(length(late))*100]; % successful early and late
      clear valid_pref mean_perf 
    end
end

mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated1 = mean_perf_consolidated;

tOutcome = cell2mat(tOutcome(:));
tOutcome1 = tOutcome;

mean_perf_consolidated = [];
tOutcome = [];

for s=bad_rec%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%intersection{s}%length(bmiBlocks)-1

      load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Performance_stats_early_late.mat']);
      mean_perf_consolidated{s,n} = mean_perf(1:2);
      early = valid_perf(1:floor(length(valid_perf)/3));
      late  = valid_perf(end-floor(length(valid_perf)/3)+1:end); 

      tOutcome{s,n} = [sum(early<15)/sum(length(early))*100,sum(late<15)/sum(length(late))*100]; % successful early and late
      clear valid_pref mean_perf 
    end
end

mean_perf_consolidated = cell2mat(mean_perf_consolidated(:));
mean_perf_consolidated2 = mean_perf_consolidated;

tOutcome = cell2mat(tOutcome(:));
tOutcome2 = tOutcome;

%% get reach data
reach_success_g = [];
count = 1; 

%all stroke rats
% for i=2:10
%     reach_success{i} =  data(days(count:count+length(intersection{stroke(i-1)})-1)+6,i);
%     days(count:count+length(intersection{stroke(i-1)})-1)
%     count = count + length(intersection{stroke(i-1)});
% end

% good recovery rats only
for i=2:7
    reach_success_g{i} =  data(days(count:count+length(intersection{good_rec(i-1)})-1)+6,i);
    days(count:count+length(intersection{good_rec(i-1)})-1)
    count = count + length(intersection{good_rec(i-1)});
end

reach_success_b = [];
count = 1; 

% bad recovery rats only
for i=2:4
    reach_success_b{i} =  data(days(count:count+length(intersection{bad_rec(i-1)})-1)+6,i);
    days(count:count+length(intersection{bad_rec(i-1)})-1)
    count = count + length(intersection{bad_rec(i-1)});
end

% days = reshape(days,1,size(days,1)*size(days,2));
reach_success_g = cell2mat(reach_success_g(:));
reach_success_b = cell2mat(reach_success_b(:));
% Correlation of reach success with BMI performance

perc1 = (mean_perf_consolidated1(:,1)-mean_perf_consolidated1(:,2))/mean_perf_consolidated1(:,1);
perc2 = (tOutcome1(:,2)-tOutcome1(:,1))/tOutcome1(:,2);

perc3 = perc1(:,5);%perc1(:,5);%perc1(:,8);
perc4 = perc2(:,1);%perc2(:,1);
% idx  = isnan(reach_success_g);
% reach_success_g(idx) = [];
% perc3(idx) = [];
% perc4(idx) = [];


perc5 = (mean_perf_consolidated2(:,1)-mean_perf_consolidated2(:,2))/mean_perf_consolidated2(:,1);
perc6 = (tOutcome2(:,2)-tOutcome2(:,1))/tOutcome2(:,2);

perc7 = perc5(:,4);%perc1(:,5);%perc1(:,8);
perc8 = perc6(:,6);%perc2(:,1);
% idx  = isnan(reach_success_b);
% reach_success_b(idx) = [];
% perc7(idx) = [];
% perc8(idx) = [];

% 
% figure();
% scatter(reach_success_g,perc3,'k','filled')
% [r,p] = corrcoef(reach_success,perc3)
% set(gca, 'TickDir', 'out')
% figure();
% scatter(reach_success_g,perc4,'k','filled')
% [r,p] = corrcoef(reach_success,perc4)
% set(gca, 'TickDir', 'out')



%% linear mixed effect model analysis 

rat_id = [];
for col=1:(size(intersection,2))
    if( ~isempty(find(col==good_rec)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end

%%
tbl = [[mean_perf_consolidated1(:,1);mean_perf_consolidated1(:,2)],[zeros(length(mean_perf_consolidated1),1);...
    ones(length(mean_perf_consolidated1),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)

%%
rat_id = [];
for col=1:(size(intersection,2))
    if( ~isempty(find(col==bad_rec)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end

%%
tbl = [[mean_perf_consolidated2(:,1);mean_perf_consolidated2(:,2)],[zeros(length(mean_perf_consolidated2),1);...
    ones(length(mean_perf_consolidated2),1)],[rat_id';rat_id']];
tbl = array2table(tbl,'VariableNames',{'time','IS','RatID'});
formula = 'time ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)


%% linear mixed model

rat_id = [];


for col=1:(size(intersection,2))
    if( ~isempty(find(col==good_rec)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end

for col=1:(size(intersection,2))
    if( ~isempty(find(col==bad_rec)) )
        rat_id = [rat_id, repelem(col,size(intersection{col},2))];
    end
end


tbl = [[mean_perf_consolidated1(:,1);mean_perf_consolidated2(:,1)],[mean_perf_consolidated1(:,2);mean_perf_consolidated2(:,2)],...
    [reach_success_g;reach_success_b],[zeros(length(mean_perf_consolidated1),1);...
    ones(length(mean_perf_consolidated2),1)],[rat_id']];
tbl = array2table(tbl,'VariableNames',{'time','success','reach','stroke_rec','RatID'});
formula = 'reach ~ time + success + stroke_rec + (time + success + stroke_rec | RatID)';
lme_R2 = fitlme(tbl,formula)

%%
figure;
plot(nanmean(data(3:5,2:10),2)); hold on;
plot(nanmean(data(7:39,2:10),2)); hold on;
% shadederror(nanmean(data(7:end,2:8),2)),nanstd)
xlim([0,33]);

%% Shaded error bar with smoothing
figure;
shadedErrorBar(data(3:5,1), smooth(nanmean(data(3:5,2:10),2)),...
     smooth(nanstd(data(3:5,2:10),0,2)/sqrt(size(data(3:5,2:10),2)-1)),'lineProps','b'); hold on;
yline(mean(data(5,2:8)));
shadedErrorBar(data(7:39,1), smooth(nanmean(data(7:39,2:10),2)),...
     smooth(nanstd(data(7:39,2:10),0,2)/sqrt(size(data(7:39,2:10),2)-1)),'lineProps','b'); hold on;
 
 %% Shaded error bar without smoothing
figure;
shadedErrorBar(data(3:5,1), (nanmean(data(3:5,2:10),2)),...
     (nanstd(data(3:5,2:10),0,2)/sqrt(size(data(3:5,2:10),2)-1)),'lineProps','b'); hold on;
yline(mean(data(5,2:10)));
shadedErrorBar(data(7:39,1), (nanmean(data(7:39,2:10),2)),...
     (nanstd(data(7:39,2:10),0,2)/sqrt(size(data(7:39,2:10),2)-1)),'lineProps','b'); hold on;


%% Scatter plot for control data - reach without BMI
% c = linspace(1,10,7);
figure;
hold all;

cmap = parula(9);
c=1:10;
% set(gca,'colororder',parula(32))

scatter(data(6:54,1),data(6:54,12:14),100,'filled'); 
colormap(gca,"winter");
%     scatter(data([2:4,6:end],1),data([2:4,6:end],i),100,'filled'); hold on;

% scatter(data(2,3:8),repelem(0.5,6),'filled'); hold on;

ylim([0,1]);
xlim([-3,40]);

yline(mean(data(6,2:10)));
set(gca, 'TickDir', 'out')