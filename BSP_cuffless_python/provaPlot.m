clc
clear all
close all

user = "Samuele";
switch user
    case "Giorgia"
        path = "/";
    case "Samuele"
        path ="/Users/samueleravazzani/LOCAL Politecnico di Milano/Local/2 anno M/BSP LAB/";
    case "Federica"
        path = "/";
end

% dati: current dataset being utilized and modified
% dati_XX_modifc: back-up
%% load data
% this database consist of a cell array of matrices, each cell is one record part.
% In each matrix each row corresponds to one signal channel: 
% 1: PPG signal, FS=125Hz;  photoplethysmograph from fingertip
% 2: ABP signal, FS=125Hz; invasive arterial blood pressure (mmHg)
% 3: ECG signal, FS=125Hz; electrocardiogram from channel II
%so part1 is struct 1x1, 3000 cells, in each of them there are matrices
%with the three signals. part1,2,3,4 are just subsets of data, probably
%because having just one struct would be computationally too heavy
part1_load= load(path + "Part_1.mat")
part2_load= load(path + "Part_2.mat")
part3_load=load(path + "Part_3.mat")
part4_load= load(path + "Part_4.mat")

dati = [part1_load.Part_1, part2_load.Part_2, part3_load.Part_3, part4_load.Part_4];
dati_00 = dati; % Back-up copy of the original
Fs = 125;
Ts = 1/Fs;

part1=part1_load.Part_1;
%% plot first signals from first part 
% PPG_1_1= part1{1, 1}(1,:);
% ABP_1_1= part1{1, 1}(2,:);
% ECG_1_1= part1{1, 1}(3,:);
% 
% figure(1)
% plot(1:size(PPG_1_1,2), PPG_1_1)
% title('PPG')
% 
% figure(2)
% plot(1:size(ABP_1_1,2), ABP_1_1)
% title('ABP')
% 
% figure(3)
% plot(1:size(ECG_1_1,2), ECG_1_1)
% title('ECG')
% 
% figure(4)
% plot(1:1000, PPG_1_1(1:1000))
% grid on
% figure(5)
% plot(1:1000, ECG_1_1(1:1000))
% grid on
% figure(6)
% plot(1:1000, ABP_1_1(1:1000))
% grid on

%%
lunghezze = [];
counter = 0;

for j = numel(dati):-1:1 %% = prod(size(cellArray)
    matrix = dati{j}; % Get the matrix from the current cell
    lunghezze = [lunghezze size(matrix,2)];
    if size(matrix,2)<37500 || any(any(isnan(matrix))) % if the length is < 5 min or there is a NaN: drop the cell
        dati(j) = [];
    else %% counts how many cells we are keeping
        counter = counter +1;
    end
    
end

disp(counter)
dati_01_removed = dati; % Back-up copy of signals after removing unaccepted signals
%% plot first 20 signals and ECG spectrum removing the mean
for i = 1:20
    figure(i); % Create a new figure for each cell
    for j = 1:4 % 3 signals + ECG spectrum
        subplot(2,2,j); % Create a subplot for each row
        if j<4
            plot(dati{i}(j,:));% Plot the row
            grid on
        else
            SGN = (dati{1}(3,:))-mean(dati{1}(3,:));
            s = fft(SGN);
            N = length(SGN);
            freq = 0:1/(N*Ts):1/Ts-1/(N*Ts);
            half = length(SGN)/2;
            plot(freq(1:half), abs(s(1:half))/N)
            grid on
        end
        switch j
            case 1
                title("PPG")
            case 2
                title("ABP")
            case 3
                title("ECG")
            case 4
                title("ECG spectrum")
        end
    end
end

%% ECG filtering: Pan Tompkins
tic
for i = 1 %% Select cell
    % Select ECG
    ECG = dati{i}(3,:);
    [ECG_peaks, ECG_peaks_indexes, ECG_delay] = pan_tompkin(ECG, Fs, 0);

    if(i==1)
        % figure(100)
        % SGN = (dati{1}(3,:))-mean(dati{1}(3,:));
        % s = fft(SGN);
        % N = length(SGN);
        % freq = 0:1/(N*Ts):1/Ts-1/(N*Ts);
        % half = length(SGN)/2;
        % plot(freq(1:half), abs(s(1:half))/N)
        % hold on
        % grid on
        % sx = fft(ECG);
        % plot(freq(1:half), abs(sx(1:half))/N)
        

    end
    % put the filter ECG back in the row of the cell
    dati{i}(3,:) = ECG;

end

tt = toc;
fprintf('Time to filter signals: %f\n', tt);
dati_02_filtered = dati; % Back-up data after filtering
%% ECG filtering
tic
for i = 1:size(dati,2) %% Select cell
    % Select ECG
    ECG = dati{i}(3,:);
    ECG = ECG-mean(ECG); % Subtract the mean

    %%%% HIGH PASS FILTER
            
    % Filter specifications
    order = 300; % Filter order
    f_c = 0.6; % Cut-off frequency in Hz
    
    % Normalized cut-off frequency (Nyquist rate)
    Wn = f_c/(Fs/2);

    % ECG = highpass(ECG,Wn,Fs);
    
    % Create the coefficients of the FIR highpass filter
    b = fir1(order, Wn, 'high');

    % Apply the filter to the signal
    ECG = filtfilt(b, 1, ECG);


    %%%% NOTCH FILTER
    % Filter specifications
    cutoff_freq = 60; % Cut-off frequency in Hz
    
    % Normalized cut-off frequency (Nyquist rate)
    Wn = cutoff_freq/(Fs/2);
    
    % Bandwidth for the notch filter
    BW = Wn/35; % common choice, rule of thumb
    
    % Create the coefficients of the FIR notch filter
    [b, a] = iirnotch(Wn, BW);
    
    % Apply the filter to the signal
    ECG = filtfilt(b, a, ECG);

    if(i==1)
        figure(100)
        SGN = (dati{1}(3,:))-mean(dati{1}(3,:));
        s = fft(SGN);
        N = length(SGN);
        freq = 0:1/(N*Ts):1/Ts-1/(N*Ts);
        half = length(SGN)/2;
        plot(freq(1:half), abs(s(1:half))/N)
        hold on
        grid on
        sx = fft(ECG);
        plot(freq(1:half), abs(sx(1:half))/N)

    end
    % put the filter ECG back in the row of the cell
    dati{i}(3,:) = ECG;

end

tt = toc;
fprintf('Time to filter signals: %f\n', tt);
dati_02_filtered = dati; % Back-up data after filtering

% PLOT
%%% plot first 20 signals and ECG spectrum FILTER
for i = 1:20
    figure(i); % Create a new figure for each cell
    for j = 1:4 % 3 signals + ECG spectrum
        subplot(2,2,j); % Create a subplot for each row
        if j<4
            plot(dati{i}(j,:));% Plot the row
            grid on
        else
            SGN = (dati{1}(3,:))-mean(dati{1}(3,:));
            s = fft(SGN);
            N = length(SGN);
            freq = 0:1/(N*Ts):1/Ts-1/(N*Ts);
            half = length(SGN)/2;
            plot(freq(1:half), abs(s(1:half))/N)
            grid on
        end
        switch j
            case 1
                title("PPG")
            case 2
                title("ABP")
            case 3
                title("ECG")
            case 4
                title("ECG spectrum")
        end
    end
end