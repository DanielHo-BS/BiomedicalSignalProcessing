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
    fs = size(data,2) / 180;
    t = linspace(0,180,size(data,2));
    f = linspace(-fs/2,fs/2,size(data,2));

    %% Separate signals from data 
    ECG = data(6,:);  % channel 6
    PPG = data(4,:);  % channel 4
    BPW = data(8,:);  % channel 8
    
    %% Band Pass Filter
    Rp = 3;
    Rs = 40;

    % High Pass Filter
    Wp = [0.3] / (fs / 2);
    Ws = [0.1] / (fs / 2);
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "high");
    tmps1 = filtfilt(B,A,ECG);
    fprintf("N and Wp of HPF = [%d, %d]\n", N, Wp)

    % Low Pass Filter
    Wp = [50] / (fs / 2);
    Ws = [55] / (fs / 2);
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "low");
    ECG_filter = filtfilt(B,A,tmps1);
    fprintf("N and Wp of LPF = [%d, %d]\n", N, Wp)
    
    
    %% Show the result
    fg = figure();
    titles = ["Original"; "HighPass Filter"; "BandPass Filter (HPF+LPF)"];
    signals = [ECG; tmps1; ECG_filter];
    for c=1:3
        subplot(3,1,c)
        plot(t, signals(c,:));
        axis([0 10 min(signals(c,:)) max(signals(c,:))])
        title(titles(c))
        xlabel("Time (sec)")
        ylabel("ECG")
    end
    saveFigure(fg,"BandPass",fname);

    %% Fourier Analysis
    fg = figure();
    signals = [ECG; ECG_filter];
    titles = ["Original"; "BandPass"];
    for c=1:2
        subplot(3,2,c)
        plot(t, signals(c,:));
        title(titles(c))
        xlabel("Time (sec)")
        ylabel("ECG")
        subplot(3,2,c+2)
        plot(f, abs(fftshift(fft(signals(c,:)))));
        xlabel("Freqence (HZ)")
        ylabel("Amplitude")
        subplot(3,2,c+4)
        plot(f, angle(fftshift(fft(signals(c,:)))));
        xlabel("Freqence (HZ)")
        ylabel("Phase")
    end
    saveFigure(fg,"Fourier",fname);  
    close all;

end

%% Save the results
function f = saveFigure(f, keyWord, fname) 
    if ~exist(append("out/",keyWord), 'dir')
       mkdir(append("out/",keyWord))
    end
    
    saveFileName = append("out/", keyWord, "/", num2str(fname));
    saveas(f,saveFileName,"jpg")

end
