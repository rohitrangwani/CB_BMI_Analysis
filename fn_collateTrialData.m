function trial_data = fn_collateTrialData(data,events,predur,postdur)
%----------------------------------------------------------------------------
% Function to collate the LFP data into channels x time x trials format
% Author@ LG from nature med paper
% INPUT -
%   -- data: LFP signals
%   -- event: Trial start events
%   -- predur: Duration before event (in seconds)
%   -- postdur: Duration after event (in seconds)
% OUTPUT - 
%   -- trail_data: A 3D array with collated LFP signals
%----------------------------------------------------------------------------

trial_data=zeros(size(data,1),predur+postdur+1,length(events));
for i = 1:size(data,1)
    for j = 1:length(events)
        if ((events(j)-predur>0)&&(events(j)+postdur<length(data)))
            trial_data(i,:,j)=data(i,events(j)-predur:events(j)+postdur);
        end
    end
end
end
