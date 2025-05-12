function fn_plotLFPs(lfp,fs,seg)
% fn_plotLFPs enables you to visuvalize short segments of LFP signal across all channels. 
% @ Aamir Abbasi. 
%
% ----------- INPUTS
% lfp = A N x M matrix where N indicate the recorded 
%       channnels and M indicate the recoded lfp amplitude in microVolts.
% fs  = The sampling frequency of lfp acquisition. 
% seg = A vector of N length in samples designating the length of lfp segment 
%        you want to visualize 
% ----------- OUTPUTS
% Just a plot!
  
h = figure; hold all;
set(0,'CurrentFigure',h);
xlim([0,5000]./fs);set(gca,'XTick',0:1:5); 
  for i=1:size(lfp,1)
    plot((1:length(lfp(i,seg)))/fs,lfp(i,seg)+(i*0.0005));
  end
end
