%% Generate binned data and plot significant units for M1

clc;clear;

%Load all session info
bmi_session_info;

before_zero = 3;
after_zero =1;
% before_zero = 1;
% after_zero =3;

tp = [];
tn = [];

res = [];
indirect = [];

for s=M1_intact%4:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%length(bmiBlocks)-1
    
        
    % Define save path
        savepath = [rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Direct_Units_Reward_AR\'];

        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*','rewards_onse*');
%         events = all_trials;
        events = rewards_onset;
          % Define save path

        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);

        valid_perf = performance(tstart{s}(n):tstop{s}(n)); 

%         indirect_unit = fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        
%         if isempty(tn_unit)
%             res{b,n} = tp_unit;
%         else
%             res{b,n} = tp_unit - tn_unit;
%         end        



        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
%         events = rewards_onset;
        
        % Get only good units
        TimeStamps = cell(size(Labels1,1),size(Labels1,2));
        for j=1:size(Labels1,1)
          for k=2:size(Labels1,2)
             if strcmp(Labels1{j,k},'good') == 1
                 TimeStamps(j,k) = TimeStamps1(j,k);
             end
          end
        end
        
        
        savepath = [rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Figs\',bmiBlocks(n).name,'\Indirect_Units_Reward_M1\'];
        [indirect_unit_M1,~,~] = fn_getPSTH_bar(TimeStamps,events,Fs,savepath,before_zero,after_zero,tstart{s}(n),tstop{s}(n));
        [indirect_base,~,~] = fn_getPSTH_bar(TimeStamps,all_trials,Fs,savepath,2,0,tstart{s}(n),tstop{s}(n));

%         indirect_unit = fn_getbinnedData(TimeStamps2,events,Fs,before_zero,after_zero,tstart{s}(n),tstop{s}(n),bin_window,0,valid_perf);
        
%         if isempty(tn_unit)
%             res{b,n} = tp_unit;
%         else
%             res{b,n} = tp_unit - tn_unit;
%         end       

        indirectM1{s,n} = indirect_unit_M1;
        base_indirect{s,n} = indirect_base;
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


%%

% indirect_ = reshape(indirect,size(indirect,1),size(indirect,2)*size(indirect,3));
indirect_M1 = reshape(indirectM1,size(indirectM1,1),size(indirectM1,2)*size(indirectM1,3));
indirect_b = reshape(base_indirect,size(base_indirect,1),size(base_indirect,2)*size(base_indirect,3));


indirectM1_cell = indirect_M1(~cellfun('isempty',indirect_M1));
indirect_cell_b = indirect_b(~cellfun('isempty',indirect_b));
% tp_cell = cell2mat(tp_cell);
% tn_cell = cell2mat(tn_cell);

indirect_data_M1 = [];
for i=1:length(indirectM1_cell)

    indirect_data_M1{i} = cell2mat(indirectM1_cell{i}(:));
    b_indirect_data{i} = cell2mat(indirect_cell_b{i}(:));
end
% indirect_data = cell2mat(indirect_data');
indirect_data_M1 = cell2mat(indirect_data_M1');
b_indirect_data = cell2mat(b_indirect_data');

%% Get significantly modulated units only

% idx1 = (4*std(indirect_data_M1(:,1:(1/2)*end),[],2) + mean(indirect_data_M1(:,1:0.5*end),2) < max(indirect_data_M1(:,.5*end:end),[],2));
% idx2 = (mean(indirect_data_M1(:,1:0.5*end),2) - 4*std(indirect_data_M1(:,1:(1/2)*end),[],2) > min(indirect_data_M1(:,.5*end:end),[],2));
% idx = idx1|idx2;
T = 800:3200;
idx1 = (4*std(b_indirect_data,[],2) + mean(b_indirect_data,2) < max(indirect_data_M1(:,T),[],2));
idx2 = (mean(b_indirect_data,2) - 4*std(b_indirect_data,[],2) > min(indirect_data_M1(:,T),[],2));
idx = idx1|idx2;
% indirect_data_M1(2.5*std(indirect_data_M1(:,1:(1/2)*end),[],2) + mean(indirect_data_M1(:,1:0.5*end),2) < max(indirect_data_M1(:,.5*end:end),[],2));

         %%
         figure;
         idx = idx1|idx2;
        indirect_data_M1_ = indirect_data_M1(idx,500:3500);
        weights = indirect_data_M1_(any(indirect_data_M1_,2),:);
        weight = normalize(weights,2,"range",[-1 1]);
        [m,idx] = max(weight(:,:),[],2);
        [ii,idx] = sort(idx); 

        weight = weight(idx,:);
        % weight = sortrows(weight,max(weight'));
        image(weight, 'CDataMapping', 'scaled')
        colorbar;
        set(gcf,'WindowState','maximized')
        set(gca, 'TickDir', 'out')
%         saveas(gcf,[rootpath,'\','IndirectUnitsM1.tiff']);
%         close;
  

%% Print all M1 robust session
robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6]};

sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};

M1 = [1,2,6,7,8];
for s=M1%length(sub)
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%length(bmiBlocks)-1
        disp(bmiBlocks(n).name)
    end
end    