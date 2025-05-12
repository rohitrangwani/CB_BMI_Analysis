function badtrials = fn_visualizeTrialData(trial_data,badchans)
%----------------------------------------------------------------------------
% Function to visualize raw LFP data and mark bad trials
% Author@ LG from nature med paper
% INPUT -
%   -- trial_data: Collated LFP signals, output from fn_collateTrialData.m
%   -- badchans: A 1D array of bad LFP channel(s) after manual inspection. 
% OUTPUT - 
%   -- badtrials: A 1D array containing the id of bad trials. 
%----------------------------------------------------------------------------

trial_data(badchans,:,:)=[]; % remove bad chans
sz=size(trial_data);
exit=0;
n=1;
badtrials=[];

fig=figure('units','normalized','outerposition',[0.1 0.1 .8 .8]);
plot(trial_data(:,:,n)');
xlim([1 sz(2)])
hold on
uicontrol('style','text','string','Bad trials','position',[50 50 80 20],'backgroundcolor',[.8 .8 .8]);
markerframes(1)=uicontrol('style','text','string','','HorizontalAlignment','left','position',[105,7,1394,40]);
hold off

while ~exit    
    title(['Trial ' num2str(n)])
    [~,~,button] = ginput(1);
    switch button
        case 32 % space
            if sum(badtrials == n)
                badtrials(badtrials==n)=[];
            else
                badtrials=[badtrials, n];                
            end
            set(markerframes(1),'string',num2str(badtrials));
        case 28 % left
            n=n-1;
        case 29 % right
            n=n+1;
        case 101 % e to exit
            exit=1;
    end
    if n>sz(3)
        exit=1;
        continue
    end
    plot(trial_data(:,:,n)');
    xlim([1 sz(2)]);
end

close(fig);