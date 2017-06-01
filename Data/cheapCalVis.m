%% Anthony Torres
% Hybrid Dev
% cheapCalVis.m
%
% Plots experimental calibration data
%
% Created: 6/1/17
% Modified: 6/1/17

%% Load data
calData = load('Data/cheapPressureTransducer.dat');
% calData = calData(1:6, :); % Remove downslope data
voltNegErr = 0.1*ones(size(calData(:,2)));
voltPosErr = 0.1*ones(size(calData(:,2)));
pressNegErr = 0.5*ones(size(calData(:,1)));
pressPosErr = 0.4*ones(size(calData(:,1)));

%% Plot data
figure;
hold on;
errorbar(calData(:,1), calData(:,2), voltNegErr, voltPosErr, pressNegErr, pressPosErr, '-o')

grid minor;
xlabel('Pressure (psi)'); ylabel('Voltage (V)');