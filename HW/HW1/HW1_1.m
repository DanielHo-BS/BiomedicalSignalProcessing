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
    
    %% Separate signals from data 
    ECG = data(6,:);  % channel 6
    PPG = data(4,:);  % channel 4
    BPW = data(8,:);  % channel 8
    
    %% Plot the signals
    signals = [BPW;PPG;ECG];
    titles = ["BPW","PPG","ECG"];
    f = figure();
    for c = 1:3
        subplot(3,1,c)
        plot(t,signals(c,:))
        title(titles(c))
        xlabel("Time(sec)")
        %axis([0 10 min(signals(c,:)) max(signals(c,:))])
    end
    
    %% Save the results
    if ~exist("out", 'dir')
       mkdir("out")
    end
    saveFileName = append("out/","HW1_180sec_",num2str(fname));
    saveas(f,saveFileName,"jpg")
    
end