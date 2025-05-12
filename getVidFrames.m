%% Getting frames for analysis

clc;clear;
tstart = {[0,0, 1, 6, 14, 5], [ 0, 0, 1, 6,  13, 5],...
    [ 0,0,0, 44, 31, 1, 15,35,8,46,13,11,38,1],...
    [ 0,0,34, 8,52, 46,21,103,20,26],...
    [ 0,0,57, 7, 8, 50,42,55,9,20]...
    [ 0,0,0,12, 9, 5, 11,14]...
    [ 0,0,12, 60, 35, 44,2,9,19]...
    [0,0,0,0,1,23,1,24,1,21]...
    [0,0, 0,9,56,10 86, 20,1]...
    [0,0, 20,58,60,1]...
    [0,0,0, 66, 12,39,2, 1, 29,14, 19, 45,17,8,  2, 61, 4, 50, 38,29, 5, 1, 15],...
    [0,0,0, 9, 6, 11, 4, 45,19, 11, 6, 1, 20,43, 17,23],...
    [0,0,0,44, 27, 57, 20, 4, 1,  33, 1 ]};

tstop  = {[0,0, 94, 106, 114, 106],[ 0, 0, 104, 105, 114, 105],...
    [ 0,0,0, 149, 91, 100, 59,138,96,139,112,33,96,93],...
    [ 0,0,133, 68, 125, 86,103,151,262,109],...
    [ 0,0,128,123,72, 154,139,111,60,50]...
    [0,0, 0,200,120,105, 134,89]...
    [0,0, 74,83,84, 80,81,103,85]...
    [0,0,0,0,116,45,43,70,32,82]...
    [0,0, 0,109,88,82,145,101,150]...
    [0,0, 80,91,103,107]...
    [0,0, 0,109,68,61,10,50,83,127,88,101,55,80,65,117,47,159,76,90,80,35,57 ],...
    [0,0,0,85,80,93,100,82,140,55,100,38,60,120,60,82],...
    [0,0,0,62, 83, 75, 48, 61,141,113,35]};

sub = {'I096','I107','I110','I111','I112', 'I122', 'I116', 'I117', 'I127', 'I128','I154','I160','I161'};

stroke= [4:5,9:13];
intact = [1:13];
intact = setdiff(intact, stroke);

Fs = 1.017252624511719e+03; 

disp('running...');
rootpath = 'Z:\Rohit\BMI_Data\';
vid_rootpath = 'Z:\Cb_BMI\';

robust_session = {[4,5],[3,4,6],[5,7,8,9,10,13,14],[3,4,5,8,9,10],[5,6,7,9],...
    [4,6,8],[5:9],[8],[7:8],[3,5,6],[8,14,15,16,18,19],[7,9,11,14],[8,9,10]};

count = 21;

for s=11:13%1:10
    s
    bmiBlocks = dir(strcat(rootpath,sub{s},'\Data\'));
%     rand_idx = randperm(length(robust_session{s}),1);
    for n=robust_session{s}
        
        %pick randonm trials 
        random_trials  = randi([tstart{s}(n),tstop{s}(n)],1,2);
        % Load vid
        root_dir = [vid_rootpath, sub{s},'\RV2_Data\', bmiBlocks(n).name];
        cd(root_dir);
        vid_file = dir('*.avi');
        vR = VideoReader(vid_file.name);

        
        
        frameRate = vR.FrameRate; % Get frame rate
        videoDuration = vR.Duration; % Get video duration
        totFrames = frameRate*videoDuration; % Calculate the total number of frames

        load([rootpath, sub{s},'\Data\', bmiBlocks(n).name,'\Events_Performance_PSTH.mat']);

        savedir = [rootpath, '\Test_vids\', bmiBlocks(n).name,'\Trial_vid\'];
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
        
        trial_vid1 = VideoWriter('vid1.avi');
        trial_vid2 = VideoWriter('vid2.avi'); 
        
        for tr=1:length(random_trials)
            clear frames vW
            frames = read(vR,[round((rewards_onset(random_trials(tr))/Fs-2.2)*vR.FrameRate),round((rewards_onset(random_trials(tr))/Fs+0.2)*vR.FrameRate)]);
            vW = VideoWriter([savedir,'Trial_',num2str(count),'.avi'],'Motion JPEG AVI');
            vW.FrameRate = vR.FrameRate;
            open(vW);
            writeVideo(vW,frames);
            close(vW);
            count = count+1;
        end
    end
end
runTime = toc(start);
disp(['done! time elapsed (minutes) - ', num2str(runTime/60)]);      
