%% HW1 Channel Separation
close all; clear; clc;
    
%% Load data 1~10 (csv file)
for fname = 1:1
    %% Read the csv file to table
    % data type: float 
    % channel number: 7
    filePath = append("../Data/LDF/",num2str(fname),".csv");
    T = readtable(filePath);
    
    % preview table
    T(1:5,:)
    % size of table
    size(T)
    
    % Calculate time and freqency of data
    % time: 0 ~ 60 sec
    t = linspace(0,60,size(T,1));
    fs = size(T,1) / 60;
    f = linspace(-fs/2,fs/2,size(T,1));
    
    % Separate signals from data 
    ECG = (T.ECG)';
    LDF_L = (T.left_LDF)';
    LDF_R = (T.right_LDF)';
    %% Show the signals
    figure;
    subplot(311),plot(t,ECG),xlabel('Time(sec)'),ylabel('ECG')
    subplot(312),plot(t,LDF_L),xlabel('Time(sec)'),ylabel('LDF_L')
    subplot(313),plot(t,LDF_R),xlabel('Time(sec)'),ylabel('LDF_R')
    
    %% find the R-Peak of ECG 
%     %% 0~4 sec / 4~12 sec / 12~60 sec
%     ECG_4 = ECG(t<4);
%     t_4 = t(t<4);
%     ECG_4_12 = ECG(t>=4 & t<12);
%     t_4_12 = t(t>=4 & t<12);
%     ECG_12 = ECG(t>=12);
%     t_12 = t(t>=12);
%     % Show signals 
%     figure;
%     subplot(411), hold on, xlabel('Time(sec)'), ylabel('ECG')
%     plot(t_4,ECG_4),plot(t_4_12,ECG_4_12),plot(t_12,ECG_12)
%     subplot(412),plot(t_4,ECG_4),xlabel('Time(sec)'),ylabel('ECG1')
%     subplot(413),plot(t_4_12,ECG_4_12),xlabel('Time(sec)'),ylabel('ECG2')
%     subplot(414),plot(t_12,ECG_12),xlabel('Time(sec)'),ylabel('ECG3')


    %% Band Pass Filter for ECG
    % High Pass Filter
%     ECG_filter1 = highPass(ECG_4, fs, 0.6, 0.3);
%     ECG_filter2 = highPass(ECG_4_12, fs, 0.6, 0.3);
%     ECG_filter3 = highPass(ECG_12, fs, 0.6, 0.3);
%     ECG_filter = highPass(ECG, fs, 0.3, 0.1);


    % Low Pass Filter
%     ECG_filter1 = lowPass(ECG_filter1, fs, 45, 55);
%     ECG_filter2 = lowPass(ECG_filter2, fs, 45, 55);
%     ECG_filter3 = lowPass(ECG_filter3, fs, 45, 55);
%     ECG_filter = [ECG_filter1,ECG_filter2,ECG_filter3];
%     ECG_filter = lowPass(ECG_filter, fs, 45, 55);


    % Show signals 
%     figure;
%     subplot(411), hold on, xlabel('Time(sec)'), ylabel('ECG')
%     plot(t_4,ECG_filter1),plot(t_4_12,ECG_filter2),plot(t_12,ECG_filter3)
%     subplot(412),plot(t_4,ECG_filter1),xlabel('Time(sec)'),ylabel('ECG1')
%     subplot(413),plot(t_4_12,ECG_filter2),xlabel('Time(sec)'),ylabel('ECG2')
%     subplot(414),plot(t_12,ECG_filter3),xlabel('Time(sec)'),ylabel('ECG3')
    ECG_filter = ECG;
                
    %% Find R Peak with window
    R=400;
    R_n=[];
    while(R < length(ECG_filter))
        index1 = R - 400;
        index2 = R + 400;
        if index1 < 1
            index1 = 1;
        end
        if index2 > length(ECG_filter)
            index2 = length(ECG_filter);
        end
        [peak,indexs] = max(ECG_filter(index1:index2));
        index = max(indexs)+index1-1;
        R_n = [R_n,index];
        R = index + 800;
    end
    %% Band Pass Filter for LDF
    % High Pass Filter
    LDF_Rf = highPass(LDF_R, fs, 0.3, 0.1);
    LDF_Lf = highPass(LDF_L, fs, 0.3, 0.1);
    %Low Pass Filter
    LDF_Rf = lowPass(LDF_Rf, fs, 15, 20);
    LDF_Lf = lowPass(LDF_Lf, fs, 15, 20);
    
    %% Find Foot of LDF with ECG's R points
    F_n_R = findFP(LDF_Rf,R_n);
    F_n_L = findFP(LDF_Lf,R_n);

    %% Show the result
    fg = figure('Position', get(0, 'Screensize'));
    subplot(311), hold on
    plot(t,ECG_filter), plot(t(R_n),ECG_filter(R_n),'O')
    title("ECG's peak points"),xlim([0,10])
    subplot(312),hold on
    plot(t,LDF_Rf),  plot(t(R_n),LDF_Rf(R_n),'o') , plot(t(F_n_R),LDF_Rf(F_n_R),'x')  
    title("LDF's peak points (Right)"),xlim([0,10])
    subplot(313), hold on
    plot(t,LDF_Lf),  plot(t(R_n),LDF_Lf(R_n),'o'), plot(t(F_n_L),LDF_Lf(F_n_L),'x')     
    title("LDF's peak points (left)"),xlim([0,10])
    saveFigure(fg, "F_Peak_10", fname, false); 
    
    %% Average wave
    ECG_m = averageWave(ECG_filter,R_n);
    LDF_Rm = averageWave(LDF_Rf,F_n_R);
    LDF_Lm = averageWave(LDF_Lf,F_n_L);
    fg = figure('Position', get(0, 'Screensize'));
    subplot(321),plot(t(1:length(ECG_m)),ECG_m), title("ECG's mean wave")
    subplot(323),plot(t(1:length(LDF_Rm)),LDF_Rm), title("LDF's mean wave (right)")
    subplot(325),plot(t(1:length(LDF_Lm)),LDF_Lm), title("LDF's mean wave (left)")
    subplot(3,2,[2,4,6]), hold on, title("Mean Wave")
    plot(t(1:length(ECG_m)),ECG_m)
    FDT = floor(mean(F_n_R - R_n));
    plot(t(1+FDT:FDT+length(LDF_Rm)),LDF_Rm)
    FDT = floor(mean(F_n_L - R_n));
    plot(t(1+FDT:FDT+length(LDF_Lm)),LDF_Lm)
    legend(["ECG", "LDF(right)", "LDF(left)"])
    saveFigure(fg, "average_wave", fname, false); 
    
    %% Analysis mean wave
    FDT = floor(mean(F_n_R - R_n));
    FR_80 = find(floor(LDF_Rm) == 80);
    FRT = max(FR_80)-min(FR_80);
    fg = figure('Position', get(0, 'Screensize'));
    hold on, title("Mean Wave")
    plot(t(1:length(ECG_m)),ECG_m)
    plot(t(1+FDT:FDT+length(LDF_Rm)),LDF_Rm)
    plot([t(min(FR_80)+FDT),t(max(FR_80)+FDT)],[LDF_Rm(min(FR_80)), LDF_Rm(max(FR_80))])
    plot([t(FDT),t(FDT)],[100,0],'black')
    legend(["ECG", "LDF(right)","FRT","FDT"])

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

%% Foot window
function F_n = findFP(sig,R_n)
    F_n = zeros(size(R_n));
    n = 1;
    for R = R_n
        % setting window size [idx1:idx2]
        idx1 = R-100;
        idx2 = R+200;
        if idx2 > length(sig)
            idx2 = length(sig);
        end
        % find the minimum of PPG from the window
        F_wave = sig(idx1:idx2);
        [~,indexs] = min(F_wave);
        % save the last index of minimum
        F_n(n) = max(indexs)+idx1-1; 
        n = n +1;
    end
end

%% High Pass Filter [ws, wp]
function sig_out = highPass(sig,fs,Wp,Ws)
    Rp = 3; Rs = 40;  % decade(dB)
    Wp = Wp / (fs / 2);
    Ws = Ws / (fs / 2);
    % using cheby1 filter
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "high");
    sig_out = filtfilt(B,A,sig);
    fprintf("N and Wp of HPF = [%d, %d]\n", N, Wp)
end

%% Low Pass Filter
function sig_out = lowPass(sig,fs, Wp, Ws)
    Rp = 3; Rs = 40;  % decade(dB)
    Wp = Wp / (fs / 2);
    Ws = Ws / (fs / 2);
    % using cheby1 filter
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "low");
    sig_out = filtfilt(B,A,sig);
    fprintf("N and Wp of LPF = [%d, %d]\n", N, Wp)
end

%% Average wave
function [waves] = averageWave(sig, sig_n)
    N_N =diff(sig_n);
    N_N_max = max(N_N);
    waves = zeros([1,N_N_max+1]);
    for i=1:length(N_N)  
        % Separate wave by foot to foot
        wave = sig(sig_n(i):sig_n(i+1));
        % Normalize wave to 0~100%
        wave_norm = (wave-min(wave))/(max(wave)-min(wave))*100;
        % Align horizontally
        wave_zero = zeros([1,N_N_max - N_N(i)]);
        % Sum all waves
        waves = waves + [wave_norm, wave_zero];
    end
    % mean the wave
    waves = waves / length(N_N);
end