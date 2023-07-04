%% HW1 Channel Separation
close all; clear; clc;
    
%% Load data 1~10 (csv file)
for fname = 2:2
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


    %% Band Pass Filter
    % High Pass Filter
%     ECG_filter1 = highPass(ECG_4, fs, 0.6, 0.3);
%     ECG_filter2 = highPass(ECG_4_12, fs, 0.6, 0.3);
%     ECG_filter3 = highPass(ECG_12, fs, 0.6, 0.3);

    % Low Pass Filter
%     ECG_filter1 = lowPass(ECG_filter1, fs, 45, 55);
%     ECG_filter2 = lowPass(ECG_filter2, fs, 45, 55);
%     ECG_filter3 = lowPass(ECG_filter3, fs, 45, 55);
%     ECG_filter = [ECG_filter1,ECG_filter2,ECG_filter3];

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
    
    %% Find Foot of LDF with ECG's R points
    F_n_R = findFP(LDF_R,R_n);
    F_n_L = findFP(LDF_L,R_n);

    %% Show the result
    fg = figure('Position', get(0, 'Screensize'));
    subplot(311), hold on
    plot(t,ECG_filter), plot(t(R_n),ECG_filter(R_n),'O')
    title("ECG's peak points"),xlim([0,10])
    subplot(312),hold on
    plot(t,LDF_R),  plot(t(R_n),LDF_R(R_n),'o') , plot(t(F_n_R),LDF_R(F_n_R),'x')  
    title("LDF's peak points (Right)"),xlim([0,10])
    subplot(313), hold on
    plot(t,LDF_L),  plot(t(R_n),LDF_L(R_n),'o'), plot(t(F_n_L),LDF_L(F_n_L),'x')     
    title("LDF's peak points (left)"),xlim([0,10])
    saveFigure(fg, "F_Peak_10", fname, false); 
    
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
