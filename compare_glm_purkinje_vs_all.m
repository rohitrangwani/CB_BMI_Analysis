%% Compare GLM models for purkinje vs others

% get all R2 for GLM modesl

bmi_session_info;

val_healthy = [];
val_stroke = [];
for s=[1,2,6,7]
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat'],'all_R2');
%     load([rootpath, sub{s}, '\GLM_model_result_purkinje.mat'],'all_R2');
    
    for j=1:10
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        
        if(~isnan(r2(4)) && r2(4))
            val_healthy_p{s,j} = r2(4);    
        end
    end
    end


    load([rootpath, sub{s}, '\GLM_M1_Cb.mat'],'all_R2');
%     load([rootpath, sub{s}, '\GLM_model_result_new.mat'],'all_R2');
    for j=1:10
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        
        if(~isnan(r2(4)) && r2(4))
            val_healthy{s,j} = r2(4);    
        end
    end
    end
end 

for s=M1_stroke
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat'],'all_R2');
%     load([rootpath, sub{s}, '\GLM_model_result_purkinje.mat'],'all_R2');
    
    for j=1:10
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(~isnan(r2(4)) && r2(4))
            val_stroke_p{s,j} = r2(4);    
        end
    end
    end



    load([rootpath, sub{s}, '\GLM_M1_Cb.mat'],'all_R2');
%     load([rootpath, sub{s}, '\GLM_model_result_new.mat'],'all_R2');
    for j=1:10
    if (sig_session_bin{s}(j))
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        
        if(~isnan(r2(4)) && r2(4))
            val_stroke{s,j} = r2(4);    
        end
    end
    end
end 

val_1 = cell2mat(val_healthy(:));

val_2 = cell2mat(val_stroke(:));


val_3 = cell2mat(val_healthy_p(:));

val_4 = cell2mat(val_stroke_p(:));

%%
val1 = val_1;
val2 = val_3;
%% plot val1 and val2 comparison

figure('Color','white'); hold all;
bar([mean(val1),mean(val2)],0.6,'w');
errorbar([1,2],[mean(val1),mean(val2)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1))]);
scatter([repelem(1,length(val1)),repelem(2,length(val2))],[val1;val2],'k');
[h,p] = ttest2(val1,val2)
set(gca,'TickDir','out');
xlim([0.5,2.5]);
% ylim([0,0.7]);
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
for col=1:(size(val_healthy_p,2))
    for row=1:(size(val_healthy_p,1))
    if(val_healthy_p{row,col})
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




%%
val1 = val_2;
val2 = val_4;
%% plot val1 and val2 comparison

figure('Color','white'); hold all;
bar([mean(val1),mean(val2)],0.6,'w');
errorbar([1,2],[mean(val1),mean(val2)],[std(val1)/sqrt((size(val2,1)-1)),std(val2)/sqrt((size(val2,1)-1))]);
scatter([repelem(1,length(val1)),repelem(2,length(val2))],[val1;val2],'k');
[h,p] = ttest2(val1,val2)
set(gca,'TickDir','out');
xlim([0.5,2.5]);
ylim([0,0.7]);
%% Mixed effect analysis

rat_id = []
for col=1:(size(val_stroke,2))
    for row=1:(size(val_stroke,1))
    if(val_stroke{row,col})
        rat_id = [rat_id,row];
    end
    end
end
rat_id
for col=1:(size(val_stroke_p,2))
    for row=1:(size(val_stroke_p,1))
    if(val_stroke_p{row,col})
        rat_id = [rat_id,row+size(val_stroke,2)];
    end
    end
end
rat_id

tbl = [[val1(:,1);val2(:,1)],[zeros(length(val1),1);...
    ones(length(val2),1)],rat_id'];
tbl = array2table(tbl,'VariableNames',{'R2','IS','RatID'});
formula = 'R2 ~ IS + (IS | RatID)';
lme_R2 = fitlme(tbl,formula)

%% Compare timescale - permutation test
%%_______________________________________________________________________%%
%
% Get highest weight time lag histogram - purkinje

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
            for b=2
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
                str_bins_p{s,j,k,b}=ii;
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
            for b=2
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
                int_bins_p{s,j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
        end
        end
    end
end 


% Get highest weight time lag histogram

bmi_session_info;

bin = [5,10,25,50,100];

str_bins = [];
int_bins = [];

for s=stroke
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_result_new.mat']);
    for j=1:size(all_model,3) 
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
        for k=1:10
            for b=2
%                 figure;
                model = all_model{b,k,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

%                 set(gcf,'WindowState','maximized')
% 
%                 saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                str_bins{s,j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
        end
        end
    end
end 

for s=intact
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_model_result_new.mat']);
    for j=1:size(all_model,3)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
        for k=1:10
            for b=2
%                 figure;
                model = all_model{b,k,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

%                 set(gcf,'WindowState','maximized')
% 
%                 saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                int_bins{s,j,k,b}=ii;
%                 tot_session = tot_session + 1;
            end    
        end
        end
    end
end
 
 %% ks test - GLM id vs GLM pd 

 str_bin = str_bins( :,:,:,2);
%  str_counts = histcounts(cell2mat(str_bin(:)),'Normalization','probability');
 int_bin = int_bins(:,:,:,2);
%  int_counts = histcounts(cell2mat(int_bin(:)),'Normalization','probability');
  str_bin_p = str_bins_p( :,:,:,2);
 int_bin_p = int_bins_p(:,:,:,2);

 
%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(str_bin(:)), cell2mat(str_bin_p(:)), 10000)

%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(int_bin(:)), cell2mat(int_bin_p(:)), 10000)
 
%  bar(int_counts); hold on;
%  bar(str_counts); hold on;


 [p, observeddifference, effectsize] = permutationTest([cell2mat(int_bin(:));cell2mat(str_bin(:))], [cell2mat(int_bin_p(:));cell2mat(str_bin_p(:))], 10000)
 
%% Get highest weight time lag histogram 

%Load all session info
bmi_session_info;
str_bins = [];
int_bins = [];
for s=M1_stroke
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_M1_Cb.mat']);
    for j=1:size(all_model,2)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
%         for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                str_bins{s,j,b}=ii;
%                 tot_session = tot_session + 1;
            end    
%         end
        end
    end
end 



for s=[1,2,6,7]
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_M1_Cb.mat']);
    for j=1:size(all_model,2)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
%         for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                int_bins{s,j,b}=ii;
%                 tot_session = tot_session + 1;
            end    
%         end
        end
    end
end 

%% Get highest weight time lag histogram - purkinje

%Load all session info
bmi_session_info;
str_bins_p = [];
int_bins_p = [];
for s=M1_stroke
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat']);
    for j=1:size(all_model,2)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
%         for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                str_bins_p{s,j,b}=ii;
%                 tot_session = tot_session + 1;
            end    
%         end
        end
    end
end 



for s=[1,2,6,7]
    s
    clear all_model all_R2 r2;
    load([rootpath, sub{s}, '\GLM_1P.mat']);
    for j=1:size(all_model,2)
        avg = mean((all_R2),2);
        r2 = avg(:,:,j);
        if(r2(2)>0)
%         for k=1:10
            for b=2:4
%                 figure;
                model = all_model{b,j};
                weights = model.Coefficients.Estimate(2:end);
                weights = reshape(abs(weights),length(weights)/((100/bin(b))+1),(100/bin(b))+1);

                weights = weights(any(weights,2),:);


                weight = normalize(weights,2,'norm',inf);

                [m,idx] = nanmax(weight,[],2);
                [ii,idx] = sort(idx); 

    %         set(gcf,'WindowState','maximized')

    %         saveas(gcf,[rootpath,'\',sub{s},'\','IndirectWeights_Sesssion_',num2str(j),'_',num2str(k),'.tiff']);
%                 close;
                int_bins_p{s,j,b}=ii;
%                 tot_session = tot_session + 1;
            end    
%         end
        end
    end
end 



 %% ks test
 str_bin = str_bins(:,:,2);
 str_counts = histcounts(cell2mat(str_bin(:)),'Normalization','probability');
 int_bin = int_bins(:,:,2);
 int_counts = histcounts(cell2mat(int_bin(:)),'Normalization','probability');
 
%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(str_bin(:)), cell2mat(int_bin(:)), 1000)
 %% ks test - GLM 1i vs GLM 1p

 str_bin = str_bins(:,:,2);
%  str_counts = histcounts(cell2mat(str_bin(:)),'Normalization','probability');
 int_bin = int_bins(:,:,2);
%  int_counts = histcounts(cell2mat(int_bin(:)),'Normalization','probability');
  str_bin_p = str_bins_p(:,:,2);
 int_bin_p = int_bins_p(:,:,2);

 
%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(str_bin(:)), cell2mat(str_bin_p(:)), 10000)

%  [h,p] = kstest2(cell2mat(int_bin(:)),cell2mat(str_bin(:)))
 [p, observeddifference, effectsize] = permutationTest(cell2mat(int_bin(:)), cell2mat(int_bin_p(:)), 10000)
 
%  bar(int_counts); hold on;
%  bar(str_counts); hold on;