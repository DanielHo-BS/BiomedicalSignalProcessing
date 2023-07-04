%% HW1 Channel Separation
close all; clear; clc;
    
%% Load data 1~10 (binary file)
result_total=zeros(10,10,2);
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
    
    %% Resampling PPG
    F_F = diff(F_n);
    F_mean = round(mean(F_F));
    PPG_resample = zeros(1,(length(F_n)-1) * (F_mean)+1);
    F_n_resample = ones(size(F_n));
    t_resample = linspace(t(F_n(1)),t(F_n(end)),size(PPG_resample,2));
    fs_resample = size(PPG_resample,2) / 180;
    f_resample = linspace(-fs/2,fs/2,size(PPG_resample,2));

    for i = 1:length(F_n)-1
        x = F_n(i):F_n(i+1);
        v = PPG_filter(x);
        xq = linspace(F_n(i),F_n(i+1),F_mean+1);
        vq = interp1(x,v,xq,'spline');
        PPG_resample((i-1)*F_mean+1:i*F_mean+1) = vq;
        F_n_resample(i+1) = (i)*F_mean+1;
    end
    
    % Plot the Resample PPG
    fg = figure('Position', get(0, 'Screensize'));
    subplot(211)
    plot(t,PPG_filter)
    hold on
    plot(t(F_n), PPG_filter(F_n), 'rO')
    title("Orignal")    

    subplot(212)
    plot(t_resample,PPG_resample)
    hold on
    plot(t_resample(F_n_resample), PPG_resample(F_n_resample), 'rO')
    title("Resample")

    saveFigure(fg,"PPG_Resample",fname,true);
    %% Fourier Analysis
    
    %% -f/2 : f/2; 0 : f/2; 0 : 20 (Hz)
    fg = figure('Position', get(0, 'Screensize'));
    signals = [PPG; PPG_filter];
    titles = ["Original"; "BandPass"];
    for c=1:2
        subplot(3,3,c)
        plot(t, signals(c,:)); xlabel("Time (sec)"); ylabel("PPG");
        title(titles(c));
        subplot(3,3,c+3)
        plot(f, abs(fftshift(fft(signals(c,:))))/length(signals(c,:))); 
        xlim([0,20]);xlabel("Freqence (HZ)"); ylabel("Amplitude");
        subplot(3,3,c+6)
        plot(f, angle(fftshift(fft(signals(c,:))))*180/pi);
        xlim([0,20]);xlabel("Freqence (HZ)"); ylabel("Phase(degree)");
    end
    subplot(3,3,3)
    plot(t_resample, PPG_resample); xlabel("Time (sec)"); ylabel("PPG");
    title("Resample");
    subplot(3,3,6)
    plot(f_resample, abs(fftshift(fft(PPG_resample)))/length(PPG_resample)); 
    xlim([0,20]);xlabel("Freqence (HZ)"); ylabel("Amplitude");
    subplot(3,3,9)
    plot(f_resample, angle(fftshift(fft(PPG_resample)))*180/pi);
    xlim([0,20]);xlabel("Freqence (HZ)"); ylabel("Phase(degree)");
    %saveFigure(fg, "Fourier_PPG", fname, true);
    %saveFigure(fg, "Fourier_PPG_half", fname, true);
    saveFigure(fg, "Fourier_PPG_20", fname, true);
    
    %% Fourier_PPG_15
    fg = figure('Position', get(0, 'Screensize'));
    subplot(3,1,1)
    plot(t_resample, PPG_resample,'r'); xlabel("Time (sec)"); ylabel("PPG");
    title("Resample");
    subplot(3,1,2)
    plot(f_resample, abs(fftshift(fft(PPG_resample)))/length(PPG_resample),'r'); 
    xlim([0,15]);xlabel("Freqence (HZ)"); ylabel("Amplitude");
    subplot(3,1,3)
    plot(f_resample, angle(fftshift(fft(PPG_resample)))*180/pi,'r');
    xlim([0,15]);xlabel("Freqence (HZ)"); ylabel("Phase(degree)");
    saveFigure(fg, "Fourier_PPG_15", fname, false);
    
    %% Periodic
    result = 0;
    for i=1:length(F_n_resample)-1
    x = PPG_resample(F_n_resample(i):F_n_resample(i+1));
    tt = t_resample(F_n_resample(i):F_n_resample(i+1));
    xdft = fft(x);
    xdft = xdft(1:floor(length(x)/2+1));
    xdft_norm = xdft/length(xdft);
    result = result + xdft_norm;
    freq = 0:fs/length(x):fs/2;
    subplot(3,1,1)
    hold on;
    plot(tt,x)
    subplot(3,1,2)
    hold on;
    stem(freq,abs(xdft_norm)); xlabel("Freq(Hz)"); ylabel("Amplitude");
    xlim([0,15]);
    end
    result = result / (length(F_n_resample)-1);
	subplot(3,1,3)
    stem(freq,abs(result)); xlabel("Freq(Hz)"); ylabel("Amplitude(Mean)");
    xlim([0,15])
    saveFigure(fg,"Fourier_PPG_Mean",fname,true);
    %%
    result_total(fname,:,1)=abs(freq(1:10));
    result_total(fname,:,2)=abs(result(1:10));
    
    
end
csvwrite('results.csv',result_total)

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
