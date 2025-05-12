function badchans = fn_checkLFPChans(lfp,seq)
%----------------------------------------------------------------------------
% Function to visualize raw LFP data and mark bad channels
% Author@ AAMIR ABBASI inspired by LG from nature med paper
% INPUT -
%   -- lfp: A 2D array of LFP signals arranged in channels x samples
%   -- seq: A 1D array containing the range of samples to plot for each LFP 
%           channel
% OUTPUT - 
%   -- badchans: A 1D array containing the id of bad LFP channels. 
%----------------------------------------------------------------------------

sz=size(lfp);
exit=0;
n=1;
badchans=[];

fig=figure('units','normalized','outerposition',[0.1 0.1 .8 .8]);
h(1) = plot(lfp(n,seq));
xlim([1 length(seq)])
hold on
uicontrol('style','text','string','Bad trials','position',[50 50 80 20],'backgroundcolor',[.8 .8 .8]);
markerframes(1)=uicontrol('style','text','string','','HorizontalAlignment','left','position',[50 25 200 20]);

while ~exit    
    title(['Channel ' num2str(n)])
    [~,~,button] = ginput(1);
    switch button
        case 32 % space
            if sum(badchans == n)
                badchans(badchans==n)=[];
            else
                badchans=[badchans, n];                
            end
            set(markerframes(1),'string',num2str(badchans));
        case 28 % left
            delete(h(n));
            n=n-1;
        case 29 % right
            n=n+1;
            if(n<=size(lfp,1))
              h(n) = plot(lfp(n,seq));
              xlim([1 length(seq)]);
            end
        case 101 % e to exit
            exit=1;
    end
    if n>sz(1)
        exit=1;
        continue
    end
end
close(fig);