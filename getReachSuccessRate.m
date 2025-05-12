% Reach Performance 


%%
sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128'};

stroke= [4:5,9:10];
rootpath = 'Z:\M1_stroke\';

% first = [4,3,4,4,3];
% last = [6,6,10,10,10];
% robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9]...
%     [4,6,8],[5:9],[8,10],[7:9],[3,5,6]};

mean_perf_consolidated = [];
success_rate = zeros(length(sub),20);
for s=stroke%4:length(sub)
    s
    blocks = dir(strcat(rootpath,sub{s}));
    for n=3:length(blocks)
      if(startsWith(blocks(n).name,"Day"))
        file = fopen([rootpath, sub{s},'\', blocks(n).name,'\Results\success_rate.txt'],'r');
        if(strlength(blocks(n).name)>4)
            d = str2num(blocks(n).name(4:5))
        else
            d = str2num(blocks(n).name(4))
        end    
        formatSpec = '%f';
        data = fscanf(file,formatSpec);
        success_rate(s,d) = data(3);
        fclose(file);
      end  
    end
end

%%
rootpath = 'Z:\Rohit\BMI_Data\Results\';

save([rootpath, '\Reach_successRate.mat'], 'success_rate');
%%
figure;
plot(success_rate(stroke,:)'); hold on;
ylim([-0.05,0.5]);