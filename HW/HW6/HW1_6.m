%% HW1 Channel Separation
close all; clear; clc;
    
%% Load data 1~10 (binary file)
for fname = 1:10
    filePath = append("../Data/",num2str(fname));
    file = fopen(filePath);
    % data type: float 
    % channel number: 8
    % time: 0 ~ 180 sec
    % sampling number: length(data)/8
    data = fread(file,[8,inf],'float');
    t = linspace(0,180,size(data,2));
    fs = size(data,2) / 180;
    f = linspace(-fs/2,fs/2,size(data,2));

    %% Separate signals from data 
    ECG = data(6,:);  % channel 6
    PPG = data(4,:);  % channel 4
    BPW = data(8,:);  % channel 8
    
    %% Load signals from result
    filePath = append("../Data/out/R_Peak/",num2str(fname),".mat");
    load(filePath)  % ECG_filter, t, R_n, R_wave
%     filePath = append("../Data/out/HRV_Resample/",num2str(fname),".mat");
%     load(filePath)  % ECG_resample, t_resample, R_n_resample
    
    %% Find the foot peak
    F_n = findFP(PPG, R_n);
    %%
    fg = figure('Position', get(0, 'Screensize'));
    subplot(211)
    plot(t, PPG)
    hold on 
    plot(t(R_n), PPG(R_n),"ro")
    xlim([0 60])
    title("PPG with R Peak of ECG")
    subplot(212)
    plot(t, PPG)
    hold on 
    plot(t(F_n), PPG(F_n),"ro")
    xlim([0 60])
    title("PPG with Foot Peak")
    saveFigure(fg, "F_Peak", fname, true);
end

%% Save the results
function fg = saveFigure(fg, keyWord, fname, closeFlage)
    if ~exist(append("out/",keyWord), 'dir')
       mkdir(append("out/",keyWord))
    end
    
    saveFileName = append("out/", keyWord, "/", num2str(fname));
    saveas(fg,saveFileName,"jpg")
    if closeFlage == true
        close;
    end
end

%% Find the R Peak
function [R_wave, R_peak, R_n] = findRP(ECG,threshold)
    % Init
    R_wave = ECG;
    ECG_min = min(ECG);
    R_wave(ECG<threshold) = ECG_min;
    R_peak = zeros(size(R_wave));
    tmp_on = 0;
    tmp_off = 0;
    % Find the Peak
    for i = 2:length(R_wave)-1
        if R_wave(i) == ECG_min && R_wave(i+1) > ECG_min
            tmp_on = i;   % Find the up pulse  
        elseif R_wave(i) == ECG_min && R_wave(i-1) > ECG_min
            tmp_off = i;  % Find the down pulse
        end
        % Find the local maximum
        if tmp_off > tmp_on  && tmp_on ~= 0 && tmp_off ~= 0
            [peak,index] = max(R_wave(tmp_on:tmp_off));
            R_peak(tmp_on+index-1) = peak;
            tmp_on = 0;
            tmp_off = 0;
        end         
    end
    R_n = find(R_peak~=0);
    R_peak = R_peak(R_n);
end

%% Adaptive
function  [R_wave, R_peak, R_n, R_R, THU] = adaptiveRP(ECG,wU,wD)
    % Init
    ECG_max = max(ECG);
    THU = ECG_max*0.8 ;
    THD = ECG_max*0.1;
    % Iteration,
    % Until Number of up line == Nummber of down line
    while(1)  
        [R_wave, R_peak, R_n_U] = findRP(ECG,THU);
        [~, ~, R_n_D] = findRP(ECG,THD);
        if length(R_n_U) ~= length(R_n_D)
            W_UD = THU - THD;
            THU = THU +wU*W_UD;
            THD = THD +wD*W_UD;
            R_n = R_n_U;
        else
            R_R = diff(R_n_U);
            fprintf("R_R:\nMin:%d\tMax:%d\tMean:%d\n",min(R_R), max(R_R), mean(R_R))
            fprintf("Threshold: %d\n",THU)
            if min(R_R) < (mean(R_R) *0.5)
                THU = THU * 1.05;
                fprintf("Refind:Threshold*1.05\n")
            elseif max(R_R) > (mean(R_R) *1.5)
                THU = THU * 0.95;
                fprintf("Refind:Threshold*0.95\n")
            else
                fprintf("Find R Wave: Done!!\n\n")
                break;
            end
            [R_wave, R_peak, R_n] = findRP(ECG,THU);
        end    
    end
end

%% Hand
function [R_wave, R_peak, R_n] = handRP(ECG,t)
    % Chose the threshold by hand
    figure();
    plot(t,ECG)          
    axis([0 10 min(ECG) max(ECG)])
    threshold = ginput(1); % Clik the point
    close;
    % Find the R Wave & Peak
    [R_wave, R_peak, R_n] = findRP(ECG,threshold(2));
end

%% Foot
function F_n = findFP(PPG,R_n)
    F_n = zeros(size(R_n));
    n = 1;
    for R = R_n
        % setting window size [idx1:idx2]
        idx1 = R;
        idx2 = R+200;
        if idx2 > length(PPG)
            idx2 = length(PPG);
        end
        % find the minimum of PPG from the window
        F_wave = PPG(idx1:idx2);
        F_min = min(F_wave);
        F = find(F_wave == F_min) + R - 1;
        % save the first index of minimum
        F_n(n) = min(F);  
        n = n +1;
    end
end
