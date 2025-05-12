%% Write tstart and tstop for all event performance files
% Rohit
clc;clear;

%Load all BMI sessions info
bmi_session_info;

%%
for s=12:length(sub)
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
    for n=robust_session{s}%first(s):last(s)%length(bmiBlocks)-1
    
        
    % Define save path
        tStart= tstart{s}(n);
        tStop= tstop{s}(n);
        save([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Valid_trials.mat'],'tStart','tStop');
        
    end
end