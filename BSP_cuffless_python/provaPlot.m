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


part1=part1_load.Part_1;
%% plot first signals from first part 
PPG_1_1= part1{1, 1}(1,:);
ABP_1_1= part1{1, 1}(2,:);
ECG_1_1= part1{1, 1}(3,:);

figure(1)
plot(1:size(PPG_1_1,2), PPG_1_1)
title('PPG')

figure(2)
plot(1:size(ABP_1_1,2), ABP_1_1)
title('ABP')

figure(3)
plot(1:size(ECG_1_1,2), ECG_1_1)
title('ECG')
%%
figure(4)
plot(1:1000, PPG_1_1(1:1000))
grid on
figure(5)
plot(1:1000, ECG_1_1(1:1000))
grid on
figure(6)
plot(1:1000, ABP_1_1(1:1000))
grid on

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

%%
%% plot first signals from first part 

% Assuming CellArray is your cell array
for i = 1:20
    figure(i); % Create a new figure for each cell
    for j = 1:4 % Assuming each matrix has 3 rows
        subplot(2,2,j); % Create a subplot for each row
        if j<4
            plot(dati{i}(j,:)); % Plot the row
        else
            Ts = 1/125;
            SGN = (dati{1}(3,:));
            s = fft(SGN);
            N = length(SGN);
            freq = 0:1/(N*Ts):1/Ts-1/(N*Ts);
            half = length(SGN)/2;
            plot(freq(1:half), abs(s(1:half))/N)
        end
    end
end



% %%
% N = 273;
% S = mean(data(2:2:N,:)+data(1:2:N-1,:)); % faccio la media del SGN: è lo stesso di prima
% W = mean(data(2:2:N,:)-data(1:2:N-1,:)); % coppie: cancello la parte comune = segnale
% SNR = std(S)/std(W)


%%
