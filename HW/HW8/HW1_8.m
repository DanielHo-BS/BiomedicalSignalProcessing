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

    %% Show the result of foot
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
        
    %% Band Pass Filter with PPG
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
    
    %% Fourier Analysis
    
    %% -f/2 : f/2
    fg = figure('Position', get(0, 'Screensize'));
    signals = [PPG; PPG_filter];
    titles = ["Original"; "BandPass"];
    for c=1:2
        subplot(3,2,c)
        plot(t, signals(c,:)); xlabel("Time (sec)"); ylabel("ECG");
        title(titles(c));
        subplot(3,2,c+2)
        plot(f, abs(fftshift(fft(signals(c,:))))/length(signals(c,:))); 
        xlim([min(f),max(f)]);xlabel("Freqence (HZ)"); ylabel("Amplitude");
        subplot(3,2,c+4)
        plot(f, angle(fftshift(fft(signals(c,:))))*180/pi);
        xlim([min(f),max(f)]);;xlabel("Freqence (HZ)"); ylabel("Phase(degree)");
    end
    saveFigure(fg, "Fourier_PPG", fname, true);
    
    %% 0 : f/2
    fg = figure('Position', get(0, 'Screensize'));
    signals = [PPG; PPG_filter];
    titles = ["Original"; "BandPass"];
    for x=1:2
    xdft = fft(signals(x,:));
    xdft = xdft(1:length(PPG)/2+1);
    xdft_norm = xdft/length(xdft);
    freq = 0:fs/length(PPG):fs/2;
    
    subplot(3,2,x)
    plot(t,PPG); xlabel("Time(sec)"); ylabel("PPG");
    title(titles(x))
    subplot(3,2,x+2)
    plot(freq,abs(xdft_norm)); xlabel("Freq(Hz)"); ylabel("Amplitude");
%     plot(freq,20*log10(abs(xdft))); xlabel("Freq(Hz)"); ylabel("Amplitude(dB)");
    xlim([min(freq),max(freq)])
    subplot(3,2,x+4)
    plot(freq,angle(xdft)*180/pi); xlabel("Freq(Hz)"); ylabel("Phase(degree)");
    xlim([min(freq),max(freq)])
    end
    saveFigure(fg, "Fourier_PPG_half", fname, true);

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
