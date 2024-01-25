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
% Create a cell array
peaks = cell(2,size(dati,2));

tic
for i = 1:size(dati,2) %% Select cell  /!\ it must be     size(dati,2)
    % Select ECG
    ECG = dati{i}(3,:);
    [ECG_filtered, ECG_peaks, ECG_peaks_indexes, ECG_delay] = pan_tompkin(ECG, Fs, 0);

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
        

        figure(102)
        plot(ECG_filtered)
        grid on
        hold on
        plot(ECG_peaks_indexes,ECG_peaks*max(ECG_filtered), 'ro')
    end

    % Put the filter ECG back in the row of the cell
    dati{i}(3,:) = ECG_filtered; % save the filtered ECG 
    
    % Save peaks
    peaks{1,i} = [ECG_peaks; ECG_peaks_indexes];
end

tt = toc;
fprintf('Time to filter ECG signals: %f\n', tt);
dati_02_filtered = dati; % Back-up data after filtering


%% PPG and ABP filtering

tic
for i = 1 %% Select cell  /!\ it must be     size(dati,2)
    % Select ECG
    PPG = dati{i}(1,:); % 1st row -> PPG
    PPG_norm = PPG/max(PPG);

    deriv1 = diff(PPG)/Ts;
    deriv1 = deriv1 /max(deriv1);
    deriv2 = diff(deriv1)/Ts;
    deriv2 = deriv2/max(deriv2);

    if(i==1)
        figure(111)
        subplot(3,1,1)
        plot(PPG_norm)
        hold on
        yline(0)
        subplot(3,1,2)
        plot(deriv1)
        hold on
        yline(0)
        subplot(3,1,3)
        plot(deriv2)
        hold on
        yline(0)
    end
    
    % Define frequency ranges
    vascular_resistance_range = [0.01, 0.15]; % LF
    
    % Apply bandpass filter for vascular range
    PPG_vascular = bandpass(PPG_fft, vascular_resistance_range);
    
    % Detect peaks using Daubechies 3 wavelet
    [PPG_maxima, PPG_max_locs] = findpeaks(PPG, 'MinPeakDistance', 100, 'MinPeakHeight', 0.5);

    % Save filtered signal

    % Save peaks
    peaks{1,i} = [ECG_peaks; ECG_peaks_indexes];
end

tt = toc;
fprintf('Time to filter signals: %f\n', tt);
dati_02_filtered = dati; % Back-up data after filtering
