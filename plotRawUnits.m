% Plot direct and indirect units examples
%


%% plot direct neurons example - Early, mid and late tr example, firing across trials, autocorrelelogram

clear; close all;
bmi_session_info;


for s=11%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(Waves_tp{i,j}))
                    

                    len = length(events);
                    figure;
                    min = 20;
                    T1 = floor(1 + (len/3-1).*rand(1,1)); 
                    T2 =  floor(len/3 + (2*len/3-len/3).*rand(1,1)); 
                    T3 =  floor(2*len/3 + (len-2*len/3).*rand(1,1)); 
%                     idx = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs;
%                     TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(idx);
%                     Waves_tp{i,j} = Waves_tp{i,j}(idx,:);
                    subplot(6,11,[1,2,12,13]);
                    if(rewards_onset(T1) == 0)
                        tr_idx  = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} < (events(T1)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T1)/Fs);
                    end

%                     tr_idx = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} <= (events(T1)/Fs +15);
                    if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                        wav = Waves_tp{i,j}(tr_idx,1:20);
                        plot([1:20],wav(1:50,1:20),'color','k'); hold on;
                    end
                    title(num2str(T1));set(gca,'TickDir','out');
                    subplot(6,11,[23,24,34,35]);
                    if(rewards_onset(T2) == 0)
                        tr_idx  = TimeStamps_tp{i,j} >= events(T2)/Fs & TimeStamps_tp{i,j} < (events(T2)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tp{i,j} >= events(T2)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T2)/Fs);
                    end
                
                   if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                        wav = Waves_tp{i,j}(tr_idx,1:20);
                        plot([1:20],wav(1:50,1:20),'color','k'); hold on;
                    end
                    title(num2str(T2));set(gca,'TickDir','out');
                    subplot(6,11,[45,46,56,57]);
                    if(rewards_onset(T3) == 0)
                        tr_idx  = TimeStamps_tp{i,j} >= events(T3)/Fs & TimeStamps_tp{i,j} < (events(T3)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tp{i,j} >= events(T3)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T3)/Fs);
                    end
                   if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                        wav = Waves_tp{i,j}(tr_idx,1:20);
                        plot([1:20],wav(1:50,1:20),'color','k'); hold on;
                    end
                    title(num2str(T3)); set(gca,'TickDir','out');
                    
                    subplot(6,11,[3:8,14:19,25:30,36:41,47:52,58:63]);
                    for tr=1:len
                        tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (events(tr)/Fs +15);
                        spikes = TimeStamps_tp{i,j}(tr_idx)-TimeStamps_tp{i,j}(find(tr_idx,1))+1;
                        scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
                    end
                    set(gca,'TickDir','out');

                    subplot(6,11,[9:11,20:22,31:33]);
                    idx = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs;
                    TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(idx);                  
                    histogram(diff(TimeStamps_tp{i,j}),'BinWidth',0.0001);
                    xlim([0,0.1]);
                    title('ISI histogram');set(gca,'TickDir','out');
                    
                    subplot(6,11,[42:44,53:55,64:66]);
                    autocorrelogram(TimeStamps_tp{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');
%                     subplot(1,4,2);
%                     plot([1:20],Waves_tp{i,j}(end/3-50:end/3,1:20),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(80:100,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,3);
%                     plot([1:20],Waves_tp{i,j}(2*end/3:2*end/3+50,1:20),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(end-20:end,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,4);
%                     hist(diff(TimeStamps_tp{i,j}(2*end/3:end,:)),1000); xlim([0,0.1]);
% %                     ylim([-5e-5,5e-5]);

                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(Waves_tn{i,j}))

                    figure;
                    T1 = floor(1 + (len/3-1).*rand(1,1)); 
                    T2 =  floor(len/3 + (2*len/3-len/3).*rand(1,1)); 
                    T3 =  floor(2*len/3 + (len-2*len/3).*rand(1,1)); 
%                     TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
%                     Waves_tn{i,j} = Waves_tn{i,j}((TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs),:);
                    subplot(6,11,[1,2,12,13]);
                    if(rewards_onset(T1) == 0)
                        tr_idx  = TimeStamps_tn{i,j} >= events(T1)/Fs & TimeStamps_tn{i,j} < (events(T1)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tn{i,j} >= events(T1)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T1)/Fs);
                    end
                    plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
                    title(num2str(T1));set(gca,'TickDir','out');
                    subplot(6,11,[23,24,34,35]);
                    if(rewards_onset(T2) == 0)
                        tr_idx  = TimeStamps_tn{i,j} >= events(T2)/Fs & TimeStamps_tn{i,j} < (events(T2)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tn{i,j} >= events(T2)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T2)/Fs);
                    end
                    plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
                    title(num2str(T2));set(gca,'TickDir','out');
                    subplot(6,11,[45,46,56,57]);
                    if(rewards_onset(T3) == 0)
                        tr_idx  = TimeStamps_tn{i,j} >= events(T3)/Fs & TimeStamps_tn{i,j} < (events(T3)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tn{i,j} >= events(T3)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T3)/Fs);
                    end
                    plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
                    title(num2str(T3));set(gca,'TickDir','out');
                    
                    
                    subplot(6,11,[3:8,14:19,25:30,36:41,47:52,58:63]);
                    for tr=1:len
                        tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (events(tr)/Fs +15);
                        spikes = TimeStamps_tn{i,j}(tr_idx)-TimeStamps_tn{i,j}(find(tr_idx,1))+1;
                        scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
                    end
                    set(gca,'TickDir','out');
                    
                    subplot(6,11,[9:11,20:22,31:33]);
                    idx = TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs;
                    TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(idx);                  
                    histogram(diff(TimeStamps_tn{i,j}),'BinWidth',0.0001);
                    xlim([0,0.1]);
                    title('ISI histogram');set(gca,'TickDir','out');


                    subplot(6,11,[42:44,53:55,64:66]);
                    autocorrelogram(TimeStamps_tn{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');
                end
            end
         end
    
    end
end

%% plot direct neurons example - Early, mid and late tr example - mean sem, firing across trials, autocorrelelogram

clear; close all;
bmi_session_info;


for s=11%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(Waves_tp{i,j}))
                    

                    len = length(events);
                    figure;

%                     idx = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs;
%                     TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(idx);
%                     Waves_tp{i,j} = Waves_tp{i,j}(idx,:);
                    subplot(6,11,[1,2,12,13]);
%                     if(rewards_onset(T1) == 0)
%                         tr_idx  = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} < (events(T1)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T1)/Fs);
%                     end
% 
% %                     tr_idx = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} <= (events(T1)/Fs +15);
%                     if(~isempty(Waves_tp{i,j}(tr_idx,:)))
%                         wav = Waves_tp{i,j}(tr_idx,1:20);
%                         plot([1:20],wav(1:50,1:20),'color','k'); hold on;
%                     end
%                     title(num2str(T1));
                    tr_idx  = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} < (events(end)/Fs +15);
                    wav =  Waves_tp{i,j}(tr_idx,1:20);
                    shadedErrorBar([1:20],wav(1:end/3,1:20),{@mean,@std});

                    set(gca,'TickDir','out');
                    subplot(6,11,[23,24,34,35]);
                    shadedErrorBar([1:20],wav(end/3:end-end/3,1:20),{@mean,@std});
                    set(gca,'TickDir','out');
                    subplot(6,11,[45,46,56,57]);
                    shadedErrorBar([1:20],wav(end-end/3:end,1:20),{@mean,@std});
                    set(gca,'TickDir','out');
                    
                    subplot(6,11,[3:8,14:19,25:30,36:41,47:52,58:63]);
                    for tr=1:len
                        tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (events(tr)/Fs +15);
                        spikes = TimeStamps_tp{i,j}(tr_idx)-TimeStamps_tp{i,j}(find(tr_idx,1))+1;
                        scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
                    end
                    set(gca,'TickDir','out');

                    subplot(6,11,[9:11,20:22,31:33]);
                    idx = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs;
                    TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(idx);                  
                    histogram(diff(TimeStamps_tp{i,j}),'BinWidth',0.0001);
                    xlim([0,0.1]);
                    title('ISI histogram');set(gca,'TickDir','out');
                    
                    subplot(6,11,[42:44,53:55,64:66]);
                    autocorrelogram(TimeStamps_tp{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');
%                     subplot(1,4,2);
%                     plot([1:20],Waves_tp{i,j}(end/3-50:end/3,1:20),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(80:100,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,3);
%                     plot([1:20],Waves_tp{i,j}(2*end/3:2*end/3+50,1:20),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(end-20:end,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,4);
%                     hist(diff(TimeStamps_tp{i,j}(2*end/3:end,:)),1000); xlim([0,0.1]);
% %                     ylim([-5e-5,5e-5]);

                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(Waves_tn{i,j}))

                    figure;
                    T1 = floor(1 + (len/3-1).*rand(1,1)); 
                    T2 =  floor(len/3 + (2*len/3-len/3).*rand(1,1)); 
                    T3 =  floor(2*len/3 + (len-2*len/3).*rand(1,1)); 
%                     TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
%                     Waves_tn{i,j} = Waves_tn{i,j}((TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs),:);
                    subplot(6,11,[1,2,12,13]);
                    if(rewards_onset(T1) == 0)
                        tr_idx  = TimeStamps_tn{i,j} >= events(T1)/Fs & TimeStamps_tn{i,j} < (events(T1)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tn{i,j} >= events(T1)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T1)/Fs);
                    end
                    plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
                    title(num2str(T1));set(gca,'TickDir','out');
                    subplot(6,11,[23,24,34,35]);
                    if(rewards_onset(T2) == 0)
                        tr_idx  = TimeStamps_tn{i,j} >= events(T2)/Fs & TimeStamps_tn{i,j} < (events(T2)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tn{i,j} >= events(T2)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T2)/Fs);
                    end
                    plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
                    title(num2str(T2));set(gca,'TickDir','out');
                    subplot(6,11,[45,46,56,57]);
                    if(rewards_onset(T3) == 0)
                        tr_idx  = TimeStamps_tn{i,j} >= events(T3)/Fs & TimeStamps_tn{i,j} < (events(T3)/Fs +15);
                    else
                        tr_idx  = TimeStamps_tn{i,j} >= events(T3)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T3)/Fs);
                    end
                    plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
                    title(num2str(T3));set(gca,'TickDir','out');
                    
                    
                    subplot(6,11,[3:8,14:19,25:30,36:41,47:52,58:63]);
                    for tr=1:len
                        tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (events(tr)/Fs +15);
                        spikes = TimeStamps_tn{i,j}(tr_idx)-TimeStamps_tn{i,j}(find(tr_idx,1))+1;
                        scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
                    end
                    set(gca,'TickDir','out');
                    
                    subplot(6,11,[9:11,20:22,31:33]);
                    idx = TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs;
                    TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(idx);                  
                    histogram(diff(TimeStamps_tn{i,j}),'BinWidth',0.0001);
                    xlim([0,0.1]);
                    title('ISI histogram');set(gca,'TickDir','out');


                    subplot(6,11,[42:44,53:55,64:66]);
                    autocorrelogram(TimeStamps_tn{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');
                end
            end
         end
    
    end
end
%% plot direct neurons example - raw example, firing across trials, isi, autocorrelelogram

clear; close all;
bmi_session_info;


for s=4:15
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(Waves_tp{i,j}))

                    len = length(events);
                    figure;
%                     T1 = floor(1 + (len/3-1).*rand(1,1)); 
%                     T2 =  floor(len/3 + (2*len/3-len/3).*rand(1,1)); 
%                     T3 =  floor(2*len/3 + (len-2*len/3).*rand(1,1)); 
% %                     idx = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs;
% %                     TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(idx);
% %                     Waves_tp{i,j} = Waves_tp{i,j}(idx,:);
%                     subplot(1,4,1);
%                     if(rewards_onset(T1) == 0)
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} < (events(T1)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T1)/Fs);
%                     end
% %                     tr_idx = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} <= (events(T1)/Fs +15);
%                     if(~isempty(Waves_tp{i,j}(tr_idx,:)))
%                         plot([1:20],Waves_tp{i,j}(tr_idx,1:20),'color','k'); hold on;
%                     end
%                     title(num2str(T1));set(gca,'TickDir','out');
%                     subplot(6,11,[23,24,34,35]);
%                     if(rewards_onset(T2) == 0)
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T2)/Fs & TimeStamps_tp{i,j} < (events(T2)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T2)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T2)/Fs);
%                     end
%                    if(~isempty(Waves_tp{i,j}(tr_idx,:)))
%                         plot([1:20],Waves_tp{i,j}(tr_idx,1:20),'color','k'); hold on;
%                     end
%                     title(num2str(T2));set(gca,'TickDir','out');
%                     subplot(6,11,[45,46,56,57]);
%                     if(rewards_onset(T3) == 0)
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T3)/Fs & TimeStamps_tp{i,j} < (events(T3)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T3)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T3)/Fs);
%                     end
%                    if(~isempty(Waves_tp{i,j}(tr_idx,:)))
%                         plot([1:20],Waves_tp{i,j}(end/3:2*end/3,1:20),'color','k'); hold on;
%                     end
%                     title(num2str(T3)); set(gca,'TickDir','out');
                    subplot(2,4,1:2);
                    if size(Waves_tp{i,j},1) > 1500
                    plot([1:20],Waves_tp{i,j}(end/3:(end/3 + 500),1:20),'color','k'); hold on;
                    end
                    subplot(2,4,[3:4,7:8]);
                    for tr=1:len
                        tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (events(tr)/Fs +15);
                        spikes = TimeStamps_tp{i,j}(tr_idx)-TimeStamps_tp{i,j}(find(tr_idx,1))+1;
                        scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
                    end
                    set(gca,'TickDir','out');

                    subplot(2,4,5);
                    idx = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs;
                    TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(idx);                  
                    histogram(diff(TimeStamps_tp{i,j}),'BinWidth',0.0001);
                    xlim([0,0.1]);
                    title('ISI histogram');set(gca,'TickDir','out');
                    
                    subplot(2,4,6);
                    autocorrelogram(TimeStamps_tp{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');

                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tp_Cb_',num2str(i),'_',num2str(j),'.fig']);
                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tp_Cb_',num2str(i),'_',num2str(j),'.tif']);
                    close
%                     subplot(1,4,2);
%                     plot([1:20],Waves_tp{i,j}(end/3-50:end/3,1:20),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(80:100,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,3);
%                     plot([1:20],Waves_tp{i,j}(2*end/3:2*end/3+50,1:20),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(end-20:end,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,4);
%                     hist(diff(TimeStamps_tp{i,j}(2*end/3:end,:)),1000); xlim([0,0.1]);
% %                     ylim([-5e-5,5e-5]);

                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(Waves_tn{i,j}))

                    figure;
%                     T1 = floor(1 + (len/3-1).*rand(1,1)); 
%                     T2 =  floor(len/3 + (2*len/3-len/3).*rand(1,1)); 
%                     T3 =  floor(2*len/3 + (len-2*len/3).*rand(1,1)); 
% %                     TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
% %                     Waves_tn{i,j} = Waves_tn{i,j}((TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs),:);
%                     subplot(6,11,[1,2,12,13]);
%                     if(rewards_onset(T1) == 0)
%                         tr_idx  = TimeStamps_tn{i,j} >= events(T1)/Fs & TimeStamps_tn{i,j} < (events(T1)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tn{i,j} >= events(T1)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T1)/Fs);
%                     end
%                     plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
%                     title(num2str(T1));set(gca,'TickDir','out');
%                     subplot(6,11,[23,24,34,35]);
%                     if(rewards_onset(T2) == 0)
%                         tr_idx  = TimeStamps_tn{i,j} >= events(T2)/Fs & TimeStamps_tn{i,j} < (events(T2)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tn{i,j} >= events(T2)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T2)/Fs);
%                     end
%                     plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
%                     title(num2str(T2));set(gca,'TickDir','out');
%                     subplot(6,11,[45,46,56,57]);
%                     if(rewards_onset(T3) == 0)
%                         tr_idx  = TimeStamps_tn{i,j} >= events(T3)/Fs & TimeStamps_tn{i,j} < (events(T3)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tn{i,j} >= events(T3)/Fs & TimeStamps_tn{i,j} < (rewards_onset(T3)/Fs);
%                     end
%                     plot([1:20],Waves_tn{i,j}(tr_idx,1:20),'color','k'); hold on;
%                     title(num2str(T3));set(gca,'TickDir','out');
%                     
%                     
%                     subplot(6,11,[3:8,14:19,25:30,36:41,47:52,58:63]);
%                     for tr=1:len
%                         tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (events(tr)/Fs +15);
%                         spikes = TimeStamps_tn{i,j}(tr_idx)-TimeStamps_tn{i,j}(find(tr_idx,1))+1;
%                         scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
%                     end
%                     set(gca,'TickDir','out');
%                     
%                     subplot(6,11,[9:11,20:22,31:33]);
%                     idx = TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs;
%                     TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(idx);                  
%                     histogram(diff(TimeStamps_tn{i,j}),'BinWidth',0.0001);
%                     xlim([0,0.1]);
%                     title('ISI histogram');set(gca,'TickDir','out');
% 
% 
%                     subplot(6,11,[42:44,53:55,64:66]);
%                     autocorrelogram(TimeStamps_tn{i,j},0.001,0.05);
%                     title('autocorrelogram');set(gca,'TickDir','out');



                    subplot(2,4,1:2);
                    if size(Waves_tn{i,j},1) > 1500
                    plot([1:20],Waves_tn{i,j}(end/3:(end/3 + 500),1:20),'color','k'); hold on;
                    end
                    subplot(2,4,[3:4,7:8]);
                    for tr=1:len
                        tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (events(tr)/Fs +15);
                        spikes = TimeStamps_tn{i,j}(tr_idx)-TimeStamps_tn{i,j}(find(tr_idx,1))+1;
                        scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
                    end
                    set(gca,'TickDir','out');

                    subplot(2,4,5);
                    idx = TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs;
                    TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(idx);                  
                    histogram(diff(TimeStamps_tn{i,j}),'BinWidth',0.0001);
                    xlim([0,0.1]);
                    title('ISI histogram');set(gca,'TickDir','out');
                    
                    subplot(2,4,6);
                    autocorrelogram(TimeStamps_tn{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');

                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tn_Cb_',num2str(i),'_',num2str(j),'.fig']);
                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tn_Cb_',num2str(i),'_',num2str(j),'.tif']);
                    close
                end
            end
         end
    
    end
end
%% plot indirect neurons example - raw example, firing across trials, isi, autocorrelelogram

clear; close all;
bmi_session_info;


for s=1%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_Cb.mat'])
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     Waves2 = Waves2';
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves2,1)
            for j=1:size(Waves2,2)
                if(~isempty(Waves2{i,j}))

                    len = length(events);
                    figure;
%                     T1 = floor(1 + (len/3-1).*rand(1,1)); 
%                     T2 =  floor(len/3 + (2*len/3-len/3).*rand(1,1)); 
%                     T3 =  floor(2*len/3 + (len-2*len/3).*rand(1,1)); 
% %                     idx = TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs;
% %                     TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(idx);
% %                     Waves_tp{i,j} = Waves_tp{i,j}(idx,:);
%                     subplot(1,4,1);
%                     if(rewards_onset(T1) == 0)
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} < (events(T1)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T1)/Fs);
%                     end
% %                     tr_idx = TimeStamps_tp{i,j} >= events(T1)/Fs & TimeStamps_tp{i,j} <= (events(T1)/Fs +15);
%                     if(~isempty(Waves_tp{i,j}(tr_idx,:)))
%                         plot([1:75],Waves_tp{i,j}(tr_idx,1:75),'color','k'); hold on;
%                     end
%                     title(num2str(T1));set(gca,'TickDir','out');
%                     subplot(6,11,[23,24,34,35]);
%                     if(rewards_onset(T2) == 0)
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T2)/Fs & TimeStamps_tp{i,j} < (events(T2)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T2)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T2)/Fs);
%                     end
%                    if(~isempty(Waves_tp{i,j}(tr_idx,:)))
%                         plot([1:75],Waves_tp{i,j}(tr_idx,1:75),'color','k'); hold on;
%                     end
%                     title(num2str(T2));set(gca,'TickDir','out');
%                     subplot(6,11,[45,46,56,57]);
%                     if(rewards_onset(T3) == 0)
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T3)/Fs & TimeStamps_tp{i,j} < (events(T3)/Fs +15);
%                     else
%                         tr_idx  = TimeStamps_tp{i,j} >= events(T3)/Fs & TimeStamps_tp{i,j} < (rewards_onset(T3)/Fs);
%                     end
%                    if(~isempty(Waves_tp{i,j}(tr_idx,:)))
%                         plot([1:75],Waves_tp{i,j}(end/3:2*end/3,1:75),'color','k'); hold on;
%                     end
%                     title(num2str(T3)); set(gca,'TickDir','out');
                    subplot(2,4,1:2);
                    plot([1:75],Waves2{i,j}(end/3:(end/3 + 500),1:75),'color','k'); hold on;
                    
                    subplot(2,4,[3:4,7:8]);
                    for tr=1:len
                        tr_idx  = TimeStamps2{i,j} >= events(tr)/Fs & TimeStamps2{i,j} < (events(tr)/Fs +15);
                        spikes = TimeStamps2{i,j}(tr_idx)-TimeStamps2{i,j}(find(tr_idx,1))+1;
                        scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
                    end
                    set(gca,'TickDir','out');

                    subplot(2,4,5);
                    idx = TimeStamps2{i,j} >= events(1)/Fs & TimeStamps2{i,j} <= events(end)/Fs;
                    TimeStamps2{i,j} = TimeStamps2{i,j}(idx);                  
                    histogram(diff(TimeStamps2{i,j}),'BinWidth',0.0001);
                    xlim([0,0.1]);
                    title('ISI histogram');set(gca,'TickDir','out');
                    
                    subplot(2,4,6);
                    autocorrelogram(TimeStamps2{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');

                    saveas(gcf,[rootpath,'Results\IndirectUnits\',bmiBlocks(n).name,'_IndirectUnit_Cb_',num2str(i),'_',num2str(j),'.fig']);
                    saveas(gcf,[rootpath,'Results\IndirectUnits\',bmiBlocks(n).name,'_IndirectUnit_Cb_',num2str(i),'_',num2str(j),'.tif']);
                    close;

                    
%                     subplot(1,4,2);
%                     plot([1:75],Waves_tp{i,j}(end/3-50:end/3,1:75),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(80:100,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,3);
%                     plot([1:20],Waves_tp{i,j}(2*end/3:2*end/3+50,1:20),'color','k'); hold on;
% %                     plot([1:30],mean(Waves_tp{i,j}(end-20:end,:),1),'color','k','linewidth',3);
%                     ylim([-5e-4,5e-4]);
%                     
%                     subplot(1,4,4);
%                     hist(diff(TimeStamps_tp{i,j}(2*end/3:end,:)),1000); xlim([0,0.1]);
% %                     ylim([-5e-5,5e-5]);

                end
            end
        end
    
    end
end

%% plot direct neurons example - raw example, firing across trials, isi, autocorrelelogram

clear; close all;
bmi_session_info;

for s=7
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(Waves_tp{i,j}))

                    len = length(events);
                    figure;

                    subplot(1,2,1);
                    if size(Waves_tp{i,j},1) > 1500
                    plot([1:20],Waves_tp{i,j}(end/3:(end/3 + 100),1:20),'color','k'); hold on;
                    end
                    
                    subplot(1,2,2);
                    autocorrelogram(TimeStamps_tp{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');

                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tp_s_Cb_',num2str(i),'_',num2str(j),'.fig']);
                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tp_s_Cb_',num2str(i),'_',num2str(j),'.tif']);
                    close

                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(Waves_tn{i,j}))

                    figure;

                    subplot(1,2,1);
                    if size(Waves_tn{i,j},1) > 1500
                    plot([1:20],Waves_tn{i,j}(end/3:(end/3 + 100),1:20),'color','k'); hold on;
                    end
                    set(gca,'TickDir','out');
                    
                    subplot(1,2,2);
                    autocorrelogram(TimeStamps_tn{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');

                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tn_s_Cb_',num2str(i),'_',num2str(j),'.fig']);
                    saveas(gcf,[rootpath,'Results\DirectUnits\',bmiBlocks(n).name,'_Tn_s_Cb_',num2str(i),'_',num2str(j),'.tif']);
                    close
                end
            end
         end
    
    end
end

%% plot indirect neurons example - raw example, firing across trials, isi, autocorrelelogram - M1

clear; close all;
bmi_session_info;


for s=12%14:15%1:15%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=14%intersection{s}(end)%first(s):last(s)%length(bmiBlocks)-1
    
     if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat'], 'file')   
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_M1.mat'])
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     Waves1 = Waves1';
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=20%4%1:size(Waves1,1)
            for j=12%4%1:size(Waves1,2)
                if(~isempty(Waves1{i,j}))

                    len = length(events);
                    figure;

                    subplot(1,2,1);
                    if size(Waves1{i,j},1) > 500
                    plot([1:40],Waves1{i,j}(end/3:(end/3 + 100),11:50),'color','k'); hold on;
                    end
%                     subplot(2,4,[3:4,7:8]);
%                     for tr=1:len
%                         tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (events(tr)/Fs +15);
%                         spikes = TimeStamps1{i,j}(tr_idx)-TimeStamps1{i,j}(find(tr_idx,1))+1;
%                         scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
%                     end
                    set(gca,'TickDir','out');
% 
%                     subplot(2,4,5);
%                     idx = TimeStamps1{i,j} >= events(1)/Fs & TimeStamps1{i,j} <= events(end)/Fs;
%                     TimeStamps2{i,j} = TimeStamps1{i,j}(idx);                  
%                     histogram(diff(TimeStamps1{i,j}),'BinWidth',0.0001);
%                     xlim([0,0.1]);
%                     title('ISI histogram');set(gca,'TickDir','out');
                    
                    subplot(1,2,2);
                    autocorrelogram(TimeStamps1{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');

                    saveas(gcf,[rootpath,'Results\IndirectUnitM1\',bmiBlocks(n).name,'_IndirectUnit_new_M1_',num2str(i),'_',num2str(j),'.fig']);
                    saveas(gcf,[rootpath,'Results\IndirectUnitM1\',bmiBlocks(n).name,'_IndirectUnit_new_M1_',num2str(i),'_',num2str(j),'.tif']);
                    close;

                    
                end
            end
        end
     end
    
    end
end

%% plot indirect neurons example - raw example, firing across trials, isi, autocorrelelogram - Cb

clear; close all;
bmi_session_info;


for s=3
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}(3)%first(s):last(s)%length(bmiBlocks)-1
    
     if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat'], 'file')   
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_Cb.mat'])
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     Waves2 = Waves2';
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1%1:size(Waves2,1)
            for j=20%1:size(Waves2,2)
                if(~isempty(Waves2{i,j}))

                    len = length(events);
                    figure;

                    subplot(1,2,1);
                    if size(Waves2{i,j},1) > 500
                    plot([1:50],Waves2{i,j}(end/3:(end/3 + 100),11:60),'color','k'); hold on;
                    end
%                     subplot(2,4,[3:4,7:8]);
%                     for tr=1:len
%                         tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (events(tr)/Fs +15);
%                         spikes = TimeStamps1{i,j}(tr_idx)-TimeStamps1{i,j}(find(tr_idx,1))+1;
%                         scatter(spikes,repelem(tr,length(spikes)),'|','k'); hold on;
%                     end
                    set(gca,'TickDir','out');
% 
%                     subplot(2,4,5);
%                     idx = TimeStamps1{i,j} >= events(1)/Fs & TimeStamps1{i,j} <= events(end)/Fs;
%                     TimeStamps2{i,j} = TimeStamps1{i,j}(idx);                  
%                     histogram(diff(TimeStamps1{i,j}),'BinWidth',0.0001);
%                     xlim([0,0.1]);
%                     title('ISI histogram');set(gca,'TickDir','out');
                    
                    subplot(1,2,2);
                    autocorrelogram(TimeStamps2{i,j},0.001,0.05);
                    title('autocorrelogram');set(gca,'TickDir','out');

                    saveas(gcf,[rootpath,'Results\IndirectUnits\',bmiBlocks(n).name,'_IndirectUnit_Cb_new_',num2str(i),'_',num2str(j),'.fig']);
                    saveas(gcf,[rootpath,'Results\IndirectUnits\',bmiBlocks(n).name,'_IndirectUnit_Cb_new_',num2str(i),'_',num2str(j),'.tif']);
                    close;

                    
                end
            end
        end
     end
    
    end
end


%% plot direct neurons 

close all;
bmi_session_info;


for s=3%:3%length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*');
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(Waves_tp{i,j}))

                    figure;
                    subplot(1,4,1);
                    Timestamps_tp{i,j} = TimeStamps_tp{i,j}(TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs);
                    Waves_tp{i,j} = Waves_tp{i,j}((TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs),:);
                   
                    hist(diff(TimeStamps_tp{i,j}(1:end/3,:)),1000); xlim([0,0.1]);
%                     plot([1:30],Waves_tp{i,j},'color',[.5 .5 .5]); hold on;
%                     plot([1:30],mean(Waves_tp{i,j}),'color','k','linewidth',3);

                    subplot(1,4,2);
                    plot([1:15],Waves_tp{i,j}(end/3-50:end/3,1:15),'color','k'); hold on;
%                     plot([1:30],mean(Waves_tp{i,j}(80:100,:),1),'color','k','linewidth',3);
                    ylim([-5e-4,5e-4]);
                    
                    subplot(1,4,3);
                    plot([1:15],Waves_tp{i,j}(2*end/3:2*end/3+50,1:15),'color','k'); hold on;
%                     plot([1:30],mean(Waves_tp{i,j}(end-20:end,:),1),'color','k','linewidth',3);
                    ylim([-5e-4,5e-4]);
                    
                    subplot(1,4,4);
                    hist(diff(TimeStamps_tp{i,j}(2*end/3:end,:)),1000); xlim([0,0.1]);
%                     ylim([-5e-5,5e-5]);

                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(Waves_tn{i,j}))

                    figure;
                    subplot(1,3,1);
                    TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
                    hist(diff(TimeStamps_tn{i,j}),1000); xlim([0,0.1]);
%                     plot([1:30],Waves_tn{i,j},'color',[.5 .5 .5]); hold on;
%                     plot([1:30],mean(Waves_tn{i,j}),'color','k','linewidth',3);
%                     ylim([-1e-4,1.5e-4]);
                    subplot(1,3,2);
                    plot([1:30],Waves_tn{i,j}(1:100,:),'color',[.5 .5 .5]); hold on;
                    plot([1:30],mean(Waves_tn{i,j}(1:100,:),1),'color','k','linewidth',3);
                    

                    subplot(1,3,3);
                    plot([1:30],Waves_tn{i,j}(end-1000:end,:),'color',[.5 .5 .5]); hold on;
                    plot([1:30],mean(Waves_tn{i,j}(end-1000:end,:),1),'color','k','linewidth',3);
                    
                    
%                     ylim([-1e-4,1.5e-4]);
                end
            end
         end
    
    end
end


%%
hist(diff(TimeStamps_tp{45,3}),10000);
xlim([0,0.1]);
%% compare early vs late ISI for all direct neurons

clear; close all;
bmi_session_info;
E_ISI = {};
L_ISI = {};
h  = {};
p = {};
k = {};

for s=1:length(sub)
    s
    count = 1;
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat'],'all_trial*');
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
        for i=1:size(TimeStamps_tp,1)
            for j=1:size(TimeStamps_tp,2)
                                    
                TimeStamps_tp{i,j} = TimeStamps_tp{i,j}(TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs);
                if(~isempty(TimeStamps_tp{i,j}))

                    
                    E_ISI{s,count} = (diff(TimeStamps_tp{i,j}(1:end/3,:)));
                    L_ISI{s,count} = (diff(TimeStamps_tp{i,j}(2*end/3:end,:)));
                    [h{s,count},p{s,count},k{s,count}] = kstest2(E_ISI{s,count},L_ISI{s,count},'Alpha',0.01);
                    count = count + 1;
    
                end
            end
        end

         for i=1:size(TimeStamps_tn,1)
            for j=1:size(TimeStamps_tn,2)
                % get spike timestamps within BMI session period
                TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
                if(~isempty(TimeStamps_tn{i,j}))


%                     hist(diff(TimeStamps_tn{i,j}),1000); xlim([0,0.1]);
                    E_ISI{s,count} = (diff(TimeStamps_tn{i,j}(1:end/3,:)));
                    L_ISI{s,count} = (diff(TimeStamps_tn{i,j}(2*end/3:end,:)));
                    [h{s,count},p{s,count},k{s,count}] = kstest2(E_ISI{s,count},L_ISI{s,count},'Alpha',0.01);
                    count = count + 1;
                end
            end
         end
    
    end
end


%% get ISI for indirect neurons - intact

clear; %close all;
bmi_session_info;
M1_ISI = {};
Cb_ISI = {};
h  = {};
p = {};
k = {};

for s=intact
    s

    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
         load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
         events = all_trials;
         events = events(tstart{s}(n):tstop{s}(n));

     if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat'], 'file')
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
        count = 1; 
        for i=1:size(TimeStamps1,1)
            for j=1:size(TimeStamps1,2)
                if(strcmp(Labels1{i,j},'good'))
                tr_count = 1;
                tr_isi = {};
                for tr=1:length(events)
                    if(rewards_onset(tr) == 0)
                        tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (events(tr)/Fs +15);
                    else
                        tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (rewards_onset(tr)/Fs);
                    end
                    
                    if(~isempty(TimeStamps1{i,j}(tr_idx)))
                    
                        tr_isi{tr_count} = (diff(TimeStamps1{i,j}(tr_idx)));
                        tr_count = tr_count + 1;
    
                    end
                end
                    M1_ISI{s,count} = cell2mat(tr_isi);
                    count = count + 1;
                end
                        
%                 TimeStamps1{i,j} = TimeStamps1{i,j}(TimeStamps1{i,j} >= events(1)/Fs & TimeStamps1{i,j} <= events(end)/Fs);

            end
        end
     end
     
     if exist([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat'], 'file')
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
         count = 1;
         for i=1:size(TimeStamps2,1)
            for j=1:size(TimeStamps2,2)
                if(strcmp(Labels2{i,j},'good'))
                % get spike timestamps within BMI session period
                   for tr=1:length(events)
                    if(rewards_onset(tr) == 0)
                        tr_idx  = TimeStamps2{i,j} >= events(tr)/Fs & TimeStamps2{i,j} < (events(tr)/Fs +15);
                    else
                        tr_idx  = TimeStamps2{i,j} >= events(tr)/Fs & TimeStamps2{i,j} < (rewards_onset(tr)/Fs);
                    end
                    if(~isempty(TimeStamps2{i,j}(tr_idx)))
             
                        Cb_ISI{s,count} = (diff(TimeStamps2{i,j}(tr_idx)));
                        count = count + 1;
                    end
                        
                   end 
                end
%                 TimeStamps2{i,j} = TimeStamps2{i,j}(TimeStamps2{i,j} >= events(1)/Fs & TimeStamps2{i,j} <= events(end)/Fs);

            end
         end
     end
    end
end


%% Analyze ISI
figure;
Cb_ISI1 = Cb_ISI(:);
Cb_ISI1(cellfun(@isempty,Cb_ISI1))= [];

Cb_ISI_3 = [];
for i=1:length(Cb_ISI1)
   Cb_ISI_3(i) = sum((Cb_ISI1{i})<0.003)/length(Cb_ISI1{i});
end

hist(Cb_ISI_3*100); %xlim([0,100]);
%%
figure;
M1_ISI1 = M1_ISI(:);
M1_ISI1(cellfun(@isempty,M1_ISI1))= [];

M1_ISI_3 = [];
for i=1:length(M1_ISI1)
   M1_ISI_3(i) = sum((M1_ISI1{i})<0.003)/length(M1_ISI1{i});
end
mean(M1_ISI_3)*100
hist(M1_ISI_3*100,100); %xlim([0,5]);


%% Analyze waveform correlation for direct neurons - with within trial template
 
close all;clear;
bmi_session_info;
Wave_corr = {};
Wave_p = {};

for s=stroke
    s
    count = 1;
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(TimeStamps_tp{i,j}))

%                     figure;
%                     Timestamps_tp{i,j} = TimeStamps_tp{i,j}(TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs);
%                     Waves_tp{i,j} = Waves_tp{i,j}((TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs),:);
                    if(~isempty(Waves_tp{i,j}))
 
                        Waves  = {};
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (rewards_onset(tr)/Fs);
                            end
                            if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves_tp{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves_tp{i,j}(tr_idx,:)';
                                [r,p] = corr(Waves{tr},Waves_temp{tr});
                                Wave_corr{s,count} = r;
                                Wave_p{s,count} = p;
                                count = count + 1;
                            end
                        end
%                         Waves = cell2mat(Waves);
                    
%                         if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(zscore(Early_wave(1:len,:))),nanmean(zscore(Late_wave(1:len,:))));
%                             Wave_corr{s,count} = [r;p];
%                         end
                    end
                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(TimeStamps_tn{i,j}))

%                     figure;
%                     TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
%                     Waves_tn{i,j} = Waves_tn{i,j}((TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs),:);
                    if(~isempty(Waves_tn{i,j}))
                        
                        Waves = {};         
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (rewards_onset(tr)/Fs);
                            end
                            if(~isempty(Waves_tn{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves_tn{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves_tn{i,j}(tr_idx,:)';
                                [r,p] = corr(Waves{tr},Waves_temp{tr});
                                Wave_corr{s,count} = r;
                                Wave_p{s,count} = p;
                                count = count + 1;
                            end
                        end
%                         Waves = cell2mat(Waves);
                        
%                        if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(Early_wave(1:len,:)),nanmean(Late_wave(1:len,:))); %[r,p] = corr(Early_wave(1:len,:),Late_wave(1:len,:));
%                             Wave_corr{s,count} = [r;p];
%                             count = count + 1;
%                        end
                    end
                end
            end
         end
    
    end
end
%% Plot wave correlation
r_val = [];
count  = 1;
figure;

for i = 1:size(Wave_corr,1)
    for j = 1:size(Wave_corr,2)
        if ~isempty(Wave_corr{i,j})
%             plot(mean(Wave_corr{i,j})); hold on;
            r_val(count) = mean(mean(Wave_corr{i,j}));
            count = count +1;
        end
    end
end

%
p_val = [];
count  = 1;
for i = 1:size(Wave_p,1)
    for j = 1:size(Wave_p,2)
        if ~isempty(Wave_p{i,j})
            p_val(count) = mean(mean(Wave_p{i,j}));
            count = count +1;
        end
    end
end

figure;
violinplot([p_val]);
figure;
violinplot(r_val);


%% Analyze waveform correlation for direct neurons - with across session template
 
close all;clear;
bmi_session_info;
Wave_corr = {};
Wave_p = {};

for s=intact
    s
    count = 1;
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(Waves_tp{i,j}))

%                     figure;
%                     Timestamps_tp{i,j} = TimeStamps_tp{i,j}(TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs);
%                     Waves_tp{i,j} = Waves_tp{i,j}((TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs),:);
                    if(~isempty(TimeStamps_tp{i,j}))
 
                        Waves  = {};
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (rewards_onset(tr)/Fs);
                            end
%                             if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves_tp{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves_tp{i,j}(tr_idx,:)';
%                                 [r,p] = corr(Waves{tr},Waves_temp{tr});
%                                 Wave_corr{s,count} = r;
%                                 Wave_p{s,count} = p;
%                                 count = count + 1;
%                             end
                        end
                        Waves_temp = cell2mat(Waves_temp);
                        
                        Template = nanmean(Waves_temp,2);
                        for tr=1:length(events)
                            if(~isempty(Waves{tr}))
                            [r,p] = corr(Waves{tr},Template);
                            Wave_corr{s,count} = r;
                            Wave_p{s,count} = p;
                            count = count + 1;
                            end
                        end
%                         Waves = cell2mat(Waves);
                    
%                         if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(zscore(Early_wave(1:len,:))),nanmean(zscore(Late_wave(1:len,:))));
%                             Wave_corr{s,count} = [r;p];
%                         end
                    end
                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(TimeStamps_tn{i,j}))

%                     figure;
%                     TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
%                     Waves_tn{i,j} = Waves_tn{i,j}((TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs),:);
                    if(~isempty(Waves_tn{i,j}))
                        
                        Waves = {};         
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (rewards_onset(tr)/Fs);
                            end
%                             if(~isempty(Waves_tn{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves_tn{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves_tn{i,j}(tr_idx,:)';
%                                 [r,p] = corr(Waves{tr},Waves_temp{tr});
%                                 Wave_corr{s,count} = r;
%                                 Wave_p{s,count} = p;
%                                 count = count + 1;
%                             end
                        end
                        Waves_temp = cell2mat(Waves_temp);
                        
                        Template = nanmean(Waves_temp,2);
                        for tr=1:length(events)
                            if(~isempty(Waves{tr}))
                            [r,p] = corr(Waves{tr},Template);
                            Wave_corr{s,count} = r;
                            Wave_p{s,count} = p;
                            count = count + 1;
                            end
                        end
%                         Waves = cell2mat(Waves);
                        
%                        if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(Early_wave(1:len,:)),nanmean(Late_wave(1:len,:))); %[r,p] = corr(Early_wave(1:len,:),Late_wave(1:len,:));
%                             Wave_corr{s,count} = [r;p];
%                             count = count + 1;
%                        end
                    end
                end
            end
         end
    
    end
end
%% Plot wave correlation
r_val = [];
count  = 1;
figure;

for i = 1:size(Wave_corr,1)
    for j = 1:size(Wave_corr,2)
        if ~isempty(Wave_corr{i,j})
%             plot(mean(Wave_corr{i,j})); hold on;
            r_val(count) = mean(mean(Wave_corr{i,j}));
            count = count +1;
        end
    end
end

%
p_val = [];
count  = 1;
for i = 1:size(Wave_p,1)
    for j = 1:size(Wave_p,2)
        if ~isempty(Wave_p{i,j})
            p_val(count) = mean(mean(Wave_p{i,j}));
            count = count +1;
        end
    end
end

figure;
violinplot([p_val]);
figure;
violinplot(r_val);


%% Plot wave correlation
r_val = [];
count  = 1;
figure;

for i = 1:size(Wave_corr,1)
    for j = 1:size(Wave_corr,2)
        if ~isempty(Wave_corr{i,j})
            plot(mean(Wave_corr{i,j}(1:15,:))); hold on;
            r_val(count) = mean(mean(Wave_corr{i,j}(1:15,:)));
            count = count +1;
        end
    end
end

%
p_val = [];
count  = 1;
for i = 1:size(Wave_corr,1)
    for j = 1:size(Wave_corr,2)
        if ~isempty(Wave_corr{i,j})
            p_val(count) = mean(mean(Wave_corr{i,j}(16:end,:)));
            count = count +1;
        end
    end
end

figure;
violinplot([p_val]);
figure;
violinplot(r_val);
%% Analyze waveform correlation for direct neurons - with across session template - early vs late
 
close all;clear;
bmi_session_info;
Wave_corr = [];
Wave_p = [];
count = 1;
for s=stroke
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
        for i=1:size(Waves_tp,1)
            for j=1:size(Waves_tp,2)
                if(~isempty(TimeStamps_tp{i,j}))

%                     figure;
%                     Timestamps_tp{i,j} = TimeStamps_tp{i,j}(TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs);
%                     Waves_tp{i,j} = Waves_tp{i,j}((TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs),:);
                    if(~isempty(Waves_tp{i,j}))
 
                        Waves  = {};
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps_tp{i,j} >= events(tr)/Fs & TimeStamps_tp{i,j} < (rewards_onset(tr)/Fs);
                            end
%                             if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves_tp{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves_tp{i,j}(tr_idx,:)';
%                                 [r,p] = corr(Waves{tr},Waves_temp{tr});
%                                 Wave_corr{s,count} = r;
%                                 Wave_p{s,count} = p;
%                                 count = count + 1;
%                             end
                        end
                        Waves_temp = cell2mat(Waves_temp);
                        
%                         Template = nanmean(Waves_temp,2);
%                         for tr=1:length(events)
%                             if(~isempty(Waves{tr}))
%                             [r,p] = corr(Waves{tr},Template);
%                             Wave_corr{s,n,tr} = r;
%                             Wave_p{s,n,tr} = p;
%                             count = count + 1;
%                             end
%                         end
                        
                        [r,p] = corr(nanmean(Waves_temp(:,1:floor(end/3)),2),nanmean(Waves_temp(:,2*end/3:end),2));
                        Wave_corr(count) = r;
                        Wave_p(count) = p;                         
                        count = count +1;
%                         Waves = cell2mat(Waves);
                    
%                         if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(zscore(Early_wave(1:len,:))),nanmean(zscore(Late_wave(1:len,:))));
%                             Wave_corr{s,count} = [r;p];
%                         end
                    end
                end
            end
        end

         for i=1:size(Waves_tn,1)
            for j=1:size(Waves_tn,2)
                if(~isempty(TimeStamps_tn{i,j}))

%                     figure;
%                     TimeStamps_tn{i,j} = TimeStamps_tn{i,j}(TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs);
%                     Waves_tn{i,j} = Waves_tn{i,j}((TimeStamps_tn{i,j} >= events(1)/Fs & TimeStamps_tn{i,j} <= events(end)/Fs),:);
                    if(~isempty(Waves_tn{i,j}))
                        
                        Waves = {};         
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps_tn{i,j} >= events(tr)/Fs & TimeStamps_tn{i,j} < (rewards_onset(tr)/Fs);
                            end
%                             if(~isempty(Waves_tn{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves_tn{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves_tn{i,j}(tr_idx,:)';
%                                 [r,p] = corr(Waves{tr},Waves_temp{tr});
%                                 Wave_corr{s,count} = r;
%                                 Wave_p{s,count} = p;
%                                 count = count + 1;
%                             end
                        end
                        Waves_temp = cell2mat(Waves_temp);
                        
%                         Template = nanmean(Waves_temp,2);
%                         for tr=1:length(events)
%                             if(~isempty(Waves{tr}))
%                             [r,p] = corr(Waves{tr},Template);
%                             Wave_corr{s,n,tr} = r;
%                             Wave_p{s,n,tr} = p;
%                             count = count + 1;
%                             end
%                         end
                        [r,p] = corr(nanmean(Waves_temp(:,1:floor(end/3)),2),nanmean(Waves_temp(:,2*end/3:end),2));
                        Wave_corr(count) = r;
                        Wave_p(count) = p;                         
                        count = count +1;
%                         Waves = cell2mat(Waves);
                        
%                        if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(Early_wave(1:len,:)),nanmean(Late_wave(1:len,:))); %[r,p] = corr(Early_wave(1:len,:),Late_wave(1:len,:));
%                             Wave_corr{s,count} = [r;p];
%                             count = count + 1;
%                        end
                    end
                end
            end
         end
    
    end
end

%%
hist(Wave_corr);

figure;
hist(Wave_p);
%% Plot wave correlation
r_val_E = [];
r_val_L = [];
count  = 1;
% figure;

for i = 1:size(Wave_corr,1)
    for j = 1:size(Wave_corr,2)
        if ~isempty(Wave_corr{i,j})
%             plot(mean(Wave_corr{i,j})); hold on;
            temp = Wave_corr(i,j,(~cellfun(@isempty,Wave_corr(i,j,:))));
            temp = squeeze(cellfun(@mean,temp));
            r_val_E(count) = mean(temp(1:end/3));
            r_val_L(count) = mean(temp(2*end/3:end));
            count = count +1;
        end
    end
end

%
p_val_E = [];
p_val_L = [];
count  = 1;
for i = 1:size(Wave_p,1)
    for j = 1:size(Wave_p,2)
        if ~isempty(Wave_p{i,j})
%             p_val(count) = mean(mean(Wave_p{i,j}));
            temp = Wave_p(i,j,(~cellfun(@isempty,Wave_p(i,j,:))));
            temp = squeeze(cellfun(@mean,temp));
            p_val_E(count) = mean(temp(1:end/3));
            p_val_L(count) = mean(temp(2*end/3:end));
            count = count +1;
        end
    end
end
[h,p] = ttest2(r_val_E,r_val_L)

[h,p] = ttest2(p_val_E,p_val_L)

% 
% figure;
% violinplot([p_val]);
% figure;
% violinplot(r_val);
%% Analyze waveform correlation for indirect neurons - with across session template
 
close all;clear;
bmi_session_info;
Wave_corr = {};
Wave_p = {};

for s=stroke
    s
    count = 1;
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_Cb.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
     Waves2 = Waves2';
     if size(Waves2,2)~=size(TimeStamps2,2)
     TimeStamps2(:,1) = [];
     end
        for i=1:size(Waves2,1)
            for j=1:size(Waves2,2)
                if(~isempty(Waves2{i,j}))

%                     figure;
%                     Timestamps_tp{i,j} = TimeStamps_tp{i,j}(TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs);
%                     Waves_tp{i,j} = Waves_tp{i,j}((TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs),:);
                    if(~isempty(TimeStamps2{i,j}))
 
                        Waves  = {};
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps2{i,j} >= events(tr)/Fs & TimeStamps2{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps2{i,j} >= events(tr)/Fs & TimeStamps2{i,j} < (rewards_onset(tr)/Fs);
                            end
%                             if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves2{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves2{i,j}(tr_idx,:)';
%                                 [r,p] = corr(Waves{tr},Waves_temp{tr});
%                                 Wave_corr{s,count} = r;
%                                 Wave_p{s,count} = p;
%                                 count = count + 1;
%                             end
                        end
                        Waves_temp = cell2mat(Waves_temp);
                        
                        Template = nanmean(Waves_temp,2);
                        for tr=1:length(events)
                            if(~isempty(Waves{tr}))
                            [r,p] = corr(Waves{tr},Template);
                            Wave_corr{s,count} = r;
                            Wave_p{s,count} = p;
                            count = count + 1;
                            end
                        end
%                         Waves = cell2mat(Waves);
                    
%                         if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(zscore(Early_wave(1:len,:))),nanmean(zscore(Late_wave(1:len,:))));
%                             Wave_corr{s,count} = [r;p];
%                         end
                    end
                end
            end
        end
    
    end
end
%% Plot wave correlation
r_val = [];
count  = 1;
figure;

for i = 1:size(Wave_corr,1)
    for j = 1:size(Wave_corr,2)
        if ~isempty(Wave_corr{i,j})
%             plot(mean(Wave_corr{i,j})); hold on;
            r_val(count) = mean(mean(Wave_corr{i,j}));
            count = count +1;
        end
    end
end

%
p_val = [];
count  = 1;
for i = 1:size(Wave_p,1)
    for j = 1:size(Wave_p,2)
        if ~isempty(Wave_p{i,j})
            p_val(count) = mean(mean(Wave_p{i,j}));
            count = count +1;
        end
    end
end

figure;
violinplot([p_val]);
figure;
violinplot(r_val);
%% Analyze waveform correlation for Cb indirect neurons - with across session template - early vs late
 
close all;clear;
bmi_session_info;
Wave_corr = [];
Wave_p = [];
count = 1;
tic;
for s=intact
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_Cb_new1.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));


          Waves2 = Waves2';
%      if size(Waves1,2)~=size(TimeStamps1,2)
%      TimeStamps2(:,1) = [];
%      Labels2(:,1) = [];
%      
        for i=1:size(Waves2,1)
            for j=1:size(Waves2,2)


                if(~isempty(Waves2{i,j}) && strcmp(Labels2{i,j},'good'))
                    if(~isempty(TimeStamps2{i,j}) && size(Waves2{i,j},1)==size(TimeStamps2{i,j},2))
                    Waves  = {};
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps2{i,j} >= events(tr)/Fs & TimeStamps2{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps2{i,j} >= events(tr)/Fs & TimeStamps2{i,j} < (rewards_onset(tr)/Fs);
                            end
                                Waves_temp{tr} = mean(Waves2{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves2{i,j}(tr_idx,:)';
                        end
                        Waves_temp = cell2mat(Waves_temp);

                        
                        [r,p] = corr(nanmean(Waves_temp(:,1:floor(end/3)),2),nanmean(Waves_temp(:,2*end/3:end),2));
                        Wave_corr(count) = r;
                        Wave_p(count) = p;                         
                        count = count +1;

                    end
                end
            end
        end    
    end
end
toc;
save([rootpath,'\Wave_corr_intact_new.mat'], 'Wave_corr','Wave_p','-v7.3');


%% Analyze waveform correlation for indirect neurons - with across session template - M1
 
close all;clear;
bmi_session_info;
Wave_corr = {};
Wave_p = {};

count = 1;
    
for s=12:15
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_M1.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));
     Waves1 = Waves1';
%      if size(Waves1,2)~=size(TimeStamps1,2)
     TimeStamps1(:,1) = [];
     Labels1(:,1) = [];
%      end
        for i=1:size(Waves1,1)
            for j=1:size(Waves1,2)
                if ~isempty(Waves1{i,j}) && strcmp(Labels1{i,j},'good')

%                     figure;
%                     Timestamps_tp{i,j} = TimeStamps_tp{i,j}(TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs);
%                     Waves_tp{i,j} = Waves_tp{i,j}((TimeStamps_tp{i,j} >= events(1)/Fs & TimeStamps_tp{i,j} <= events(end)/Fs),:);
                    if(~isempty(TimeStamps1{i,j}) && size(Waves1{i,j},1)==size(TimeStamps1{i,j},2))
 
                        Waves  = {};
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (rewards_onset(tr)/Fs);
                            end
%                             if(~isempty(Waves_tp{i,j}(tr_idx,:)))
                                Waves_temp{tr} = mean(Waves1{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves1{i,j}(tr_idx,:)';
%                                 [r,p] = corr(Waves{tr},Waves_temp{tr});
%                                 Wave_corr{s,count} = r;
%                                 Wave_p{s,count} = p;
%                                 count = count + 1;
%                             end
                        end
                        Waves_temp = cell2mat(Waves_temp);
                        
                        Template = nanmean(Waves_temp,2);
                        for tr=1:length(events)
                            if(~isempty(Waves{tr}))
                            [r,p] = corr(Waves{tr},Template);
                            Wave_corr{s,count} = r;
                            Wave_p{s,count} = p;
                            count = count + 1;
                            end
                        end
%                         Waves = cell2mat(Waves);
                    
%                         if(~isempty(Waves))
%                             Early_wave = Waves(1:15,1:floor(end/3));
%                             Late_wave = Waves(1:15,floor(2*end/3):end);
%                             len =  min(size(Early_wave,2),size(Late_wave,2));
%                             [r,p] = corr(Early_wave(:,1:len),Late_wave(:,1:len));
%     %                         [r,p] = corr(nanmean(zscore(Early_wave(1:len,:))),nanmean(zscore(Late_wave(1:len,:))));
%                             Wave_corr{s,count} = [r;p];
%                         end
                    end
                end
            end
        end
    
    end
end

%% Plot wave correlation
r_val = [];
count  = 1;
figure;

for i = 1:size(Wave_corr,1)
    for j = 1:size(Wave_corr,2)
        if ~isempty(Wave_corr{i,j})
%             plot(mean(Wave_corr{i,j})); hold on;
            r_val(count) = mean(mean(Wave_corr{i,j}));
            count = count +1;
        end
    end
end

%
p_val = [];
count  = 1;
for i = 1:size(Wave_p,1)
    for j = 1:size(Wave_p,2)
        if ~isempty(Wave_p{i,j})
            p_val(count) = mean(mean(Wave_p{i,j}));
            count = count +1;
        end
    end
end

figure;
violinplot([p_val]);
figure;
violinplot(r_val);

%% Analyze waveform correlation for M1 indirect neurons - with across session template - early vs late
 
close all;clear;
bmi_session_info;
Wave_corr = [];
Wave_p = [];
count = 1;
for s=M1_stroke
    s

    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_M1.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
     events = all_trials;
     events = events(tstart{s}(n):tstop{s}(n));
     rewards_onset = rewards_onset(tstart{s}(n):tstop{s}(n));


          Waves1 = Waves1';
%      if size(Waves1,2)~=size(TimeStamps1,2)
     TimeStamps1(:,1) = [];
     Labels1(:,1) = [];
%      
        for i=1:size(Waves1,1)
            for j=1:size(Waves1,2)


                if(~isempty(Waves1{i,j}) && strcmp(Labels1{i,j},'good'))
                    if(~isempty(TimeStamps1{i,j}) && size(Waves1{i,j},1)==size(TimeStamps1{i,j},2))

                    Waves  = {};
                        Waves_temp = {};
                        for tr=1:length(events)
                            if(rewards_onset(tr) == 0)
                                tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (events(tr)/Fs +15);
                            else
                                tr_idx  = TimeStamps1{i,j} >= events(tr)/Fs & TimeStamps1{i,j} < (rewards_onset(tr)/Fs);
                            end
                                Waves_temp{tr} = mean(Waves1{i,j}(tr_idx,:),1)';
                                Waves{tr} = Waves1{i,j}(tr_idx,:)';
                        end
                        Waves_temp = cell2mat(Waves_temp);

                        
                        [r,p] = corr(nanmean(Waves_temp(:,1:floor(end/3)),2),nanmean(Waves_temp(:,2*end/3:end),2));
                        Wave_corr(count) = r;
                        Wave_p(count) = p;                         
                        count = count +1;

                    end
                end
            end
        end    
    end
end

save([rootpath,'\M1_wave_corr_stroke_new.mat'], 'Wave_corr','Wave_p','-v7.3');

%% Load indirect neuron example - Cb

load('\\smb.researchcampus.csmc.edu\gulatitlab-data\Rohit\BMI_Data\I096\Data\I096-211105-114444\Timestamps_M1.mat')
%%
figure;
plot([1:30],Waves_tp{45,2}(1:100,:),'color',[.5 .5 .5]); hold on;
plot([1:30],mean(Waves_tp{45,2}(1:100,:),1),'color','k','linewidth',3);

subp
plot([1:30],Waves_tp{45,2}(end-100:end,:),'color',[.5 .5 .5]); hold on;
plot([1:30],mean(Waves_tp{45,2}(end-100:end,:),1),'color','k','linewidth',3);


%% Load indirect neuron example - M1

load('\\smb.researchcampus.csmc.edu\gulatitlab-data\Rohit\BMI_Data\I096\Data\I096-211105-114444\Timestamps_Direct.mat')
%%
figure;
plot([1:30],Waves_tp{45,2}(1:100,:),'color',[.5 .5 .5]); hold on;
plot([1:30],mean(Waves_tp{45,2}(1:100,:),1),'color','k','linewidth',3);

subp
plot([1:30],Waves_tp{45,2}(end-100:end,:),'color',[.5 .5 .5]); hold on;
plot([1:30],mean(Waves_tp{45,2}(end-100:end,:),1),'color','k','linewidth',3);


%%

%% Compute L-Ratio and Isolation Distance 

close all;clear;
bmi_session_info;
isolationDistance = [];
L_ratio = [];
tic;
count = 1;
for s=stroke
    s

    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Cb.mat'],'Labels2');
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_Cb_new1.mat']);
        Waves2 = Waves2';
        if(~isempty(Waves2))
            isolation_qual{s,n} =  compute_unit_isolation(Waves2,Labels2,2);
        end
%         for i=1:size(Waves2,1)
%             for j=1:size(Waves2,2)
%                 if (~isempty(Waves2{i,j}) && ~isempty(TimeStamps2{i,j}))
%                      spikeWaveforms = Waves2{i,j};
%                      unitIdx = TimeStamps2{i,j};
%                     
% %                     % Assume we have PCA-based spike features (rows: spikes, cols: features)
% %                     numSpikes = length(spikeTimes);
% %                     numFeatures = 3; % Number of PCA features
% %                     spikeFeatures = randn(numSpikes, numFeatures); % Simulated feature matrix
% %                     unitIdx = randi([0 1], numSpikes, 1); % 1 = unit spikes, 0 = noise spikes
% %                     
% %                     % Compute Mahalanobis distances of unit spikes to noise distribution
% %                     mu_noise = mean(spikeFeatures(unitIdx == 0, :), 1);
% %                     cov_noise = cov(spikeFeatures(unitIdx == 0, :));
% %                     mahal_dist = mahal(spikeFeatures(unitIdx == 1, :), mu_noise, cov_noise);
% %                     
% %                     % Compute Isolation Distance
% %                     n_noise = sum(unitIdx == 0);
% %                     sortedDist = sort(mahal_dist);
% %                     if n_noise < length(sortedDist)
% %                         isolationDistance(count) = sortedDist(n_noise);
% %                     else
% %                         isolationDistance(count) = Inf;
% %                     end
% %                     
% %                     % Compute L-ratio
% %                     df = numFeatures; % Degrees of freedom
% %                     p_values = 1 - chi2cdf(mahal_dist, df);
% %                     L_ratio(count) = sum(p_values) / length(unitIdx == 1);
% 
% 
%                     % Compute L-Ratio and Isolation Distance (Simulated Feature Data)
% %                     numSpikes = length(spikeTimes);
%                     % Perform PCA on Spike Waveforms
%                     numPCs = 3; % Number of principal components to use
%                     [coeff, score, ~] = pca(spikeWaveforms); % Perform PCA
%                     spikeFeatures = score(:,1:numPCs); % Extract first 3 PCA components
% %                     unitIdx = randi([0 1], numSpikes, 1); % 1 = unit spikes, 0 = noise spikes
%                     
%                     % Compute Mahalanobis distances of unit spikes to noise distribution
%                     if sum(unitIdx == 0) > numPCs % Ensure enough noise points
%                         mahal_dist = mahal(spikeFeatures(unitIdx == 1, :), spikeFeatures(unitIdx == 0, :));
%                     else
%                         mahal_dist = NaN(size(spikeFeatures(unitIdx == 1, :), 1), 1);
%                     end
%                     
%                     % Compute Isolation Distance
%                     n_noise = sum(unitIdx == 0);
%                     sortedDist = sort(mahal_dist(~isnan(mahal_dist)));
%                     if n_noise < length(sortedDist)
%                         isolationDistance(count) = sortedDist(n_noise);
%                     else
%                         isolationDistance(count) = Inf;
%                     end
%                     
%                     % Compute L-ratio
%                     df = numPCs; % Degrees of freedom
%                     p_values = 1 - chi2cdf(mahal_dist, df);
%                     L_ratio(count) = sum(p_values, 'omitnan') / sum(unitIdx == 1);
%                     count = count +1;
%                 end
%             end
%         end
        toc
    end
end
%
% save([rootpath,'\Stroke_isolation.mat'], 'isolation_qual','-v7.3');
save([rootpath,'\Stroke_isolation_new1.mat'], 'isolation_qual','-v7.3');

%%
iso_dist = [];
iso_l = [];
for i = 1:size(isolation_qual,1)
    for j = 1:size(isolation_qual,2)
        temp =  isolation_qual{i,j};
        if (length(temp)>1)
            iso_dist = [iso_dist [temp.isolation_distance]];
            iso_l = [iso_l [temp.L_ratio]];
        end
    end
end


%% Compute L-Ratio and Isolation Distance - M1

close all;clear;
bmi_session_info;
isolationDistance = [];
L_ratio = [];
tic;
count = 1;
for s=M1_stroke
    s

    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1

        
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_M1.mat'],'Labels1');
        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Waves_M1.mat']);
        Waves1 = Waves1';
        Labels1(:,1) = [];
        isolation_qual{s,n} =  compute_unit_isolation(Waves1,Labels1,2);
%         for i=1:size(Waves2,1)
%             for j=1:size(Waves2,2)
%                 if (~isempty(Waves2{i,j}) && ~isempty(TimeStamps2{i,j}))
%                      spikeWaveforms = Waves2{i,j};
%                      unitIdx = TimeStamps2{i,j};
%                     
% %                     % Assume we have PCA-based spike features (rows: spikes, cols: features)
% %                     numSpikes = length(spikeTimes);
% %                     numFeatures = 3; % Number of PCA features
% %                     spikeFeatures = randn(numSpikes, numFeatures); % Simulated feature matrix
% %                     unitIdx = randi([0 1], numSpikes, 1); % 1 = unit spikes, 0 = noise spikes
% %                     
% %                     % Compute Mahalanobis distances of unit spikes to noise distribution
% %                     mu_noise = mean(spikeFeatures(unitIdx == 0, :), 1);
% %                     cov_noise = cov(spikeFeatures(unitIdx == 0, :));
% %                     mahal_dist = mahal(spikeFeatures(unitIdx == 1, :), mu_noise, cov_noise);
% %                     
% %                     % Compute Isolation Distance
% %                     n_noise = sum(unitIdx == 0);
% %                     sortedDist = sort(mahal_dist);
% %                     if n_noise < length(sortedDist)
% %                         isolationDistance(count) = sortedDist(n_noise);
% %                     else
% %                         isolationDistance(count) = Inf;
% %                     end
% %                     
% %                     % Compute L-ratio
% %                     df = numFeatures; % Degrees of freedom
% %                     p_values = 1 - chi2cdf(mahal_dist, df);
% %                     L_ratio(count) = sum(p_values) / length(unitIdx == 1);
% 
% 
%                     % Compute L-Ratio and Isolation Distance (Simulated Feature Data)
% %                     numSpikes = length(spikeTimes);
%                     % Perform PCA on Spike Waveforms
%                     numPCs = 3; % Number of principal components to use
%                     [coeff, score, ~] = pca(spikeWaveforms); % Perform PCA
%                     spikeFeatures = score(:,1:numPCs); % Extract first 3 PCA components
% %                     unitIdx = randi([0 1], numSpikes, 1); % 1 = unit spikes, 0 = noise spikes
%                     
%                     % Compute Mahalanobis distances of unit spikes to noise distribution
%                     if sum(unitIdx == 0) > numPCs % Ensure enough noise points
%                         mahal_dist = mahal(spikeFeatures(unitIdx == 1, :), spikeFeatures(unitIdx == 0, :));
%                     else
%                         mahal_dist = NaN(size(spikeFeatures(unitIdx == 1, :), 1), 1);
%                     end
%                     
%                     % Compute Isolation Distance
%                     n_noise = sum(unitIdx == 0);
%                     sortedDist = sort(mahal_dist(~isnan(mahal_dist)));
%                     if n_noise < length(sortedDist)
%                         isolationDistance(count) = sortedDist(n_noise);
%                     else
%                         isolationDistance(count) = Inf;
%                     end
%                     
%                     % Compute L-ratio
%                     df = numPCs; % Degrees of freedom
%                     p_values = 1 - chi2cdf(mahal_dist, df);
%                     L_ratio(count) = sum(p_values, 'omitnan') / sum(unitIdx == 1);
%                     count = count +1;
%                 end
%             end
%         end
        toc
    end
end
%
save([rootpath,'\M1_Stroke_isolation_new.mat'], 'isolation_qual','-v7.3');


%% Compute L-Ratio and Isolation Distance - direct units
close all;clear;
bmi_session_info;
isolationDistance = [];
L_ratio = [];
tic;
count = 1;
for s=stroke
    s

    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=intersection{s}%first(s):last(s)%length(bmiBlocks)-1
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Timestamps_Direct.mat']);
     load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);
        
        label = ~cellfun(@isempty, Waves_tp);
     
        if sum(sum(label)) > 0
        isolation_qual_tp{s,n} =  compute_unit_isolation_direct(Waves_tp,label,2);
        end
        label = ~cellfun(@isempty, Waves_tn);
        if sum(sum(label)) > 0
        isolation_qual_tn{s,n} =  compute_unit_isolation_direct(Waves_tn,label,2);
        end
        toc
    end
end
%
% save([rootpath,'\M1_Stroke_isolation_new.mat'], 'isolation_qual','-v7.3');