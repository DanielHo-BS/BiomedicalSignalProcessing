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
    %% Band Pass Filter
    Rp = 3;
    Rs = 40;

    % High Pass Filter
    Wp = 0.3 / (fs / 2);
    Ws = 0.1 / (fs / 2);
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "high");
    tmps1 = filtfilt(B,A,PPG);
    %fprintf("N and Wp of HPF = [%d, %d]\n", N, Wp)

    % Low Pass Filter
    Wp = 50 / (fs / 2);
    Ws = 55 / (fs / 2);
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "low");
    PPG_filter = filtfilt(B,A,tmps1);
    %fprintf("N and Wp of LPF = [%d, %d]\n", N, Wp)

    %% Find the foot peak
    % force-find foots by Adaptive
    [F_wave, F_peak, F_n, F_F, THU] = adaptiveFP(PPG_filter,-0.1,0.1);
    % find foots again by window
    F_n = findFP_byR(PPG, F_n);
    F_F = diff(F_n);
    fprintf("F_F:\nMin:%d\tMax:%d\tMean:%d\n",min(F_F), max(F_F), mean(F_F))
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

%% Find the Foot
function [F_wave, F_peak, F_n] = findFP(PPG,threshold)
    % Init
    F_wave = PPG;
    F_wave(PPG>threshold) = threshold;
    F_peak = zeros(size(F_wave));
    tmp_on = 0;
    tmp_off = 0;
    % Find the foot
    for i = 2:length(F_wave)-1
        if F_wave(i) == threshold && F_wave(i+1) < threshold
            tmp_on = i;   % Find the down pulse  
        elseif F_wave(i) == threshold && F_wave(i-1) < threshold
            tmp_off = i;  % Find the up pulse
        end
        % Find the local maximum
        if tmp_off > tmp_on  && tmp_on ~= 0 && tmp_off ~= 0
            [peak,indexs] = min(F_wave(tmp_on:tmp_off));
            index = max(indexs);  % find the last one
            F_peak(tmp_on+index-1) = peak;
            tmp_on = 0;
            tmp_off = 0;
        end         
    end
    F_n = find(F_peak~=0);
end

%% Adaptive
function  [F_wave, F_peak, F_n, F_F, THU] = adaptiveFP(PPG,wU,wD)
    % Init
    PPG_min = min(PPG);
    THD = PPG_min*0.8 ;
    THU = PPG_min*0.1;
    % Iteration,
    % Until Number of up line == Nummber of down line
    while(1)  
        [F_wave, F_peak, F_n_U] = findFP(PPG,THU);
        [~, ~, F_n_D] = findFP(PPG,THD);
        if length(F_n_U) ~= length(F_n_D)
            W_UD = THU - THD;
            THU = THU +wU*W_UD;
            THD = THD +wD*W_UD;
            F_n = F_n_U;
        else
            F_F = diff(F_n_U);
            % find the errors: 
            % Too long (Foot to Foot)
            errors = find(F_F > mean(F_F)*1.5);
            while(not(isempty(errors)))
                errors = find(F_F > mean(F_F)*1.5);
                for error = errors
                    num = round(F_F(error)/mean(F_F));
                    for i = 1:num
                        mid = round(F_n(error)+ mean(F_F)*i);
                        index1 = mid-50;
                        index2 = mid+50;
                        if index1 < 1
                            index1 = 1;
                        end
                        if index2 > length(PPG)
                            index2 = length(PPG);
                        end
                        [peak,indexs] = min(PPG(index1:index2));
                        index = max(indexs)+mid-100-1;  % find the last one
                        F_n = [F_n, index];
                    end
                end
                F_n = sort(F_n);
                % Too short (Foot to Foot)
                F_F = diff(F_n);
                F_n = [F_n(F_F > mean(F_F)*0.5), F_n(end)];
                F_F = diff(F_n);
            end
            fprintf("Threshold: %d\n",THU)
            break
        end    
    end
end

%% Foot
function F_n = findFP_byR(PPG,R_n)
    F_n = zeros(size(R_n));
    n = 1;
    R_R = diff(R_n);
    num = round(R_n(1)/mean(R_R));
	for i = 1:num
        mid = round(1 + mean(R_R)*i);
        index1 = mid-50;
        index2 = mid+50;
        if index1 < 1
            index1 = 1;
        end
        if index2 > length(PPG)
            index2 = length(PPG);
        end
        [~,indexs] = min(PPG(index1:index2));
        index = max(indexs)+mid-100-1;  % find the last one
        R_n = [R_n, index];
    end
    R_n = sort(R_n);
    for R = R_n
        % setting window size [idx1:idx2]
        idx1 = R-400;
        idx2 = R+100;
        if idx1 < 1
            idx1 = 1;
        end
        if idx2 > length(PPG)
            idx2 = length(PPG);
        end
        
        % find the minimum of PPG from the window
        F_wave = PPG(idx1:idx2);
        [F_min ,F] = min(F_wave);
        F = F + idx1 - 1;
        % save the first index of minimum
        F_n(n) = max(F);  
        n = n +1;
    end
    % Too short (Foot to Foot)
    F_F = diff(F_n);
    F_n = [F_n(F_F > mean(F_F)*0.5), F_n(end)];
end
