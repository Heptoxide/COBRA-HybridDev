%% Anthony Torres
% COBRA
% chamberFoS.m
% 
% Calculates the temperature dependence of SS304, and outputs the resulting
% factor of safey

%% Initialize data
% Pulled from: http://www.aksteel.com/pdf/markets_products/stainless/austenitic/304_304l_data_bulletin.pdf
temps = [25; 204; 316; 427; 538; 649; 760; 871]; % deg C
ss304_YS = [35; 23; 20; 17; 14; 13; 11]; % ksi

% Plot yield strength vs temp