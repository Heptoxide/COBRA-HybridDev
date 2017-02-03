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

%% General motor properties
burnTime = 10; % s


%% Initialize fuel properties
a_HTPB = 0.198; % mm/s
a_Paraffin = 0.1146; % mm/s

n_HTPB = 0.325;
n_Paraffin = 0.5036;

rho_HTPB = (899.5+930)/2; % kg/m^3
rho_Paraffin = (900+930)/2; % kg/m^3

%%%%% UPDATE THESE VALUES %%%%%
mfDot = 0.038; % kg/s
moDot = 0.263; % kg/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initPortRad = linspace(convlength(0.25, 'in', 'm'), convlength(1.25, 'in', 'm'));


%% Fuel length
fuelLen_Paraffin = calcFuelLength(mfDot, moDot, rho_Paraffin, a_Paraffin, n_Paraffin, initPortRad, 1);
fuelLen_HTPB = calcFuelLength(mfDot, moDot, rho_HTPB, a_HTPB, n_HTPB, initPortRad, 1);

% Plot results
% Doing calculations in metric, but showing results in imperial because
% it's more intuitive for the US folks (who will be reviewing this)
figure;
hold on;
plot(convlength(initPortRad, 'm', 'in'), convlength(fuelLen_Paraffin, 'm', 'in'), 'b');
plot(convlength(initPortRad, 'm', 'in'), convlength(fuelLen_HTPB, 'm', 'in'), 'r');
hold off;

grid on; grid minor;
legend('Paraffin', 'HTPB');

xlabel('Initial Port Radius (in)'); ylabel('Required Fuel Length (in)');


%% Fuel mass
fuelMass_Paraffin = calcTotalFuel(moDot, rho_Paraffin, a_Paraffin, ...
                      n_Paraffin, initPortRad, 1, min(fuelLen_Paraffin), burnTime);
                  
fuelMass_HTPB = calcTotalFuel(moDot, rho_HTPB, a_HTPB, ...
                   n_HTPB, initPortRad, 1, min(fuelLen_HTPB), burnTime);

% Plot results
% Doing calculations in metric, but showing results in imperial because
% it's more intuitive for the US folks (who will be reviewing this)
figure;
hold on;
plot(convlength(initPortRad, 'm', 'in'), fuelMass_Paraffin, 'b');
plot(convlength(initPortRad, 'm', 'in'), fuelMass_HTPB, 'r');
hold off;

grid on; grid minor;
legend('Paraffin', 'HTPB', 'Location', 'SouthEast');

xlabel('Initial Port Radius (in)'); ylabel('Required Fuel Mass (kg)');


%% Outer radius
outerRadius_Paraffin = calcInstFuelRadius(moDot, a_Paraffin, n_Paraffin, initPortRad, 1, burnTime);
outerRadius_HTPB = calcInstFuelRadius(moDot, a_HTPB, n_HTPB, initPortRad, 1, burnTime);

% Plot results
% Doing calculations in metric, but showing results in imperial because
% it's more intuitive for the US folks (who will be reviewing this)
figure;
hold on;
plot(convlength(initPortRad, 'm', 'in'), convlength(outerRadius_Paraffin, 'm', 'in'), 'b');
plot(convlength(initPortRad, 'm', 'in'), convlength(outerRadius_HTPB, 'm', 'in'), 'r');
hold off;

grid on; grid minor;
legend('Paraffin', 'HTPB', 'Location', 'SouthEast');

xlabel('Initial Port Radius (in)'); ylabel('Outer Fuel Grain Radius (in)');