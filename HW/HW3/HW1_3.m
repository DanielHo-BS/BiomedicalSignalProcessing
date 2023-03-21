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
    Wp = 0.6 / (fs / 2);
    Ws = 0.3 / (fs / 2);
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "high");
    tmps1 = filtfilt(B,A,ECG);
    %fprintf("N and Wp of HPF = [%d, %d]\n", N, Wp)

    % Low Pass Filter
    Wp = 50 / (fs / 2);
    Ws = 55 / (fs / 2);
    [N, Wp]=cheb1ord(Wp, Ws, Rp, Rs);
    [B,A] = cheby1(N, Rp, Wp, "low");
    ECG_filter = filtfilt(B,A,tmps1);
    %fprintf("N and Wp of LPF = [%d, %d]\n", N, Wp)
    
    %% Find R Wave
    byHand = false;
    if byHand == true  % Hand
    [R_wave, R_peak, R_n] = handRP(ECG_filter,t);

    else  % Adaptive
        wU = -0.01;  % weight of up line
        wD = 0.01;   % weight of down line
        [R_wave, R_peak, R_n] = adaptiveRP(ECG_filter,wU,wD);
    end
    %fprintf("Num of R: %d\n",length(R_n))
    
    % Show the result
    fg = figure();
    subplot(2,1,1)
    plot(t,ECG_filter)          
    hold on
    plot(t(R_n),R_peak,'O')
    axis([0 20 min(ECG_filter) max(ECG_filter)])
    title(["Num of R: ",num2str(length(R_n))])
    
    subplot(2,1,2)
    plot(t,R_wave)
    hold on
    plot(t(R_n),R_peak,'O')
    axis([0 20 min(ECG_filter) max(ECG_filter)])
    
    saveFigure(fg,"R Peak",fname,true);
    %% Show the result
%     fg = figure();
%     titles = ["Original"; "HighPass Filter"; "BandPass Filter (HPF+LPF)"];
%     signals = [ECG; tmps1; ECG_filter];
%     for c=1:3
%         subplot(3,1,c)
%         plot(t, signals(c,:));
%         axis([0 10 min(signals(c,:)) max(signals(c,:))])
%         title(titles(c))
%         xlabel("Time (sec)")
%         ylabel("ECG")
%     end
%     saveFigure(fg,"BandPass",fname);

    %% Fourier Analysis
%     fg = figure();
%     signals = [ECG; ECG_filter];
%     titles = ["Original"; "BandPass"];
%     for c=1:2
%         subplot(3,2,c)
%         plot(t, signals(c,:));
%         title(titles(c))
%         xlabel("Time (sec)")
%         ylabel("ECG")
%         subplot(3,2,c+2)
%         plot(f, abs(fftshift(fft(signals(c,:)))));
%         xlabel("Freqence (HZ)")
%         ylabel("Amplitude")
%         subplot(3,2,c+4)
%         plot(f, angle(fftshift(fft(signals(c,:)))));
%         xlabel("Freqence (HZ)")
%         ylabel("Phase")
%     end
%     %saveFigure(fg,"Fourier",fname,true);  

end

%% Save the results
function f = saveFigure(f, keyWord, fname,closeFlage)
    if ~exist(append("out/",keyWord), 'dir')
       mkdir(append("out/",keyWord))
    end
    
    saveFileName = append("out/", keyWord, "/", num2str(fname));
    saveas(f,saveFileName,"jpg")
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
function  [R_wave, R_peak, R_n_U] = adaptiveRP(ECG,wU,wD)
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
        else
            break;
        end      
    end
    disp(THU)
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