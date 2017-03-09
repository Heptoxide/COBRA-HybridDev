%% Anthony Torres
% COBRA - Hybrid Dev
% fuel_Sens.m
%
% Created: 2/3/17
% Modified: 2/3/17
%
% To Do:
%   Add in other fuels, especially aluminized paraffin
%
% Purpose: Calculates various fuel properties for sensitivity analysis
%       length
%       mass
%       outer radius
%       thrust
%
% Limitations:
%   Only works for circular port geometries
%
% Inputs:
%   None
%
% Outputs:
%   
%

%% Clear things
clc; close all;

%% General motor properties
burnTime = 5; % s


%% Initialize fuel properties
% ABS regression constants: http://trace.tennessee.edu/cgi/viewcontent.cgi?article=3962&context=utk_gradthes

a_HTPB = 0.198; % mm/s
a_Paraffin = 0.1146; % mm/s
a_ABS = 0.417; % mm/s

n_HTPB = 0.325;
n_Paraffin = 0.5036;
n_ABS = 0.347; % mm/s

rho_HTPB = (899.5+930)/2; % kg/m^3
rho_Paraffin = (900+930)/2; % kg/m^3

% ABS densities: http://www.stelray.com/reference-tables.html
rho_ABS_70F = 1040; % kg/m^3
rho_ABS_Melt = 970; % kg/m^3

%%%%% UPDATE THESE VALUES %%%%%
mfDot = 0.038; % kg/s
moDot = 0.263; % kg/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initPortRad = linspace(convlength(0.125, 'in', 'm'), convlength(1.25, 'in', 'm'));


%% Fuel length
fuelLen_Paraffin = calcFuelLength(mfDot, moDot, rho_Paraffin, a_Paraffin, n_Paraffin, initPortRad, 1);
fuelLen_HTPB = calcFuelLength(mfDot, moDot, rho_HTPB, a_HTPB, n_HTPB, initPortRad, 1);
fuelLen_ABS_70F = calcFuelLength(mfDot, moDot, rho_ABS_70F, a_ABS, n_ABS, initPortRad, 1);
fuelLen_ABS_Melt = calcFuelLength(mfDot, moDot, rho_ABS_Melt, a_ABS, n_ABS, initPortRad, 1);

% Plot results
% Doing calculations in metric, but showing results in imperial because
% it's more intuitive for the US folks (who will be reviewing this)
figure;
hold on;
plot(convlength(initPortRad, 'm', 'in'), convlength(fuelLen_Paraffin, 'm', 'in'), 'b');
plot(convlength(initPortRad, 'm', 'in'), convlength(fuelLen_HTPB, 'm', 'in'), 'r');
plot(convlength(initPortRad, 'm', 'in'), convlength(fuelLen_ABS_70F, 'm', 'in'), '-.m');
plot(convlength(initPortRad, 'm', 'in'), convlength(fuelLen_ABS_Melt, 'm', 'in'), ':m');
vline(0.25, '--r');
hold off;

grid on; grid minor;
legend('Paraffin', 'HTPB', 'ABS 70^\circ F', 'ABS Melt');

xlabel('Initial Port Radius (in)'); ylabel('Required Fuel Length (in)');


%% Fuel mass
fuelMass_Paraffin = calcTotalFuel(moDot, rho_Paraffin, a_Paraffin, ...
                      n_Paraffin, initPortRad, 1, min(fuelLen_Paraffin), burnTime);
                  
fuelMass_HTPB = calcTotalFuel(moDot, rho_HTPB, a_HTPB, ...
                   n_HTPB, initPortRad, 1, min(fuelLen_HTPB), burnTime);

fuelMass_ABS_70F = calcTotalFuel(moDot, rho_ABS_70F, a_ABS, ...
                   n_ABS, initPortRad, 1, min(fuelLen_ABS_70F), burnTime);

fuelMass_ABS_Melt = calcTotalFuel(moDot, rho_ABS_Melt, a_ABS, ...
                   n_ABS, initPortRad, 1, min(fuelLen_ABS_Melt), burnTime);
               
% Plot results
% Doing calculations in metric, but showing results in imperial because
% it's more intuitive for the US folks (who will be reviewing this)
figure;
hold on;
plot(convlength(initPortRad, 'm', 'in'), fuelMass_Paraffin, 'b');
plot(convlength(initPortRad, 'm', 'in'), fuelMass_HTPB, 'r');
plot(convlength(initPortRad, 'm', 'in'), fuelMass_ABS_70F, '-.m');
plot(convlength(initPortRad, 'm', 'in'), fuelMass_ABS_Melt, ':m');
vline(0.25, '--r');
hold off;

grid on; grid minor;
legend('Paraffin', 'HTPB', 'ABS 70^\circ F', 'ABS Melt', 'Location', 'SouthEast');

xlabel('Initial Port Radius (in)'); ylabel('Required Fuel Mass (kg)');


%% Outer radius
outerRadius_Paraffin = calcInstFuelRadius(moDot, a_Paraffin, n_Paraffin, initPortRad, 1, burnTime);
outerRadius_HTPB = calcInstFuelRadius(moDot, a_HTPB, n_HTPB, initPortRad, 1, burnTime);
outerRadius_ABS_70F = calcInstFuelRadius(moDot, a_ABS, n_ABS, initPortRad, 1, burnTime);
outerRadius_ABS_Melt = calcInstFuelRadius(moDot, a_ABS, n_ABS, initPortRad, 1, burnTime);

% Plot results
% Doing calculations in metric, but showing results in imperial because
% it's more intuitive for the US folks (who will be reviewing this)
figure;
hold on;
plot(convlength(initPortRad, 'm', 'in'), convlength(outerRadius_Paraffin, 'm', 'in'), 'b');
plot(convlength(initPortRad, 'm', 'in'), convlength(outerRadius_HTPB, 'm', 'in'), 'r');
plot(convlength(initPortRad, 'm', 'in'), convlength(outerRadius_ABS_70F, 'm', 'in'), '-.m');
plot(convlength(initPortRad, 'm', 'in'), convlength(outerRadius_ABS_Melt, 'm', 'in'), ':m');
vline(0.25, '--r');
hold off;

grid on; grid minor;
legend('Paraffin', 'HTPB', 'ABS 70^\circ F', 'ABS Melt', 'Location', 'SouthEast');

xlabel('Initial Port Radius (in)'); ylabel('Outer Fuel Grain Radius (in)');