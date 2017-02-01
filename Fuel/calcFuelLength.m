%% Anthony Torres
% COBRA - Hybrid Dev
% calcFuelLength.m
%
% Created: 2/1/17
% Modified: 2/1/17
%
% Purpose: Calculates the required fuel grain length with the given
% parameters. Rearrangement of RPE (15-12)
%
% Limitations:
%   Only works for circular port geometries
%
% Inputs:
%   mfDot --- initial mass flow rate of fuel (kg/s)
%
%   moDot --- mass flow rate of oxidizer (kg/s)
%
%   rhoF --- density of fuel (kg/m^3)
%         900 to 930 kg/m^3 for Paraffin
%         ~899.5 to 930 kg/m^3 for HTPB R45
%
%   a --- regression coefficient for O-F combination (mm/s)
%         0.198 mm/s for HTPB
%         0.1146 mm/s for Paraffin
%
%   n --- regression exponent coefficient for O-F combination
%         0.325 for HTPB
%         0.5036 for Paraffin
%
%   Ri --- initial port radius (m)
%
%   N --- number of circular ports
%
% Outputs:
%   L --- length of fuel grain (m)
%

function[L] = calcFuelLength(mfDot, moDot, rhoF, a, n, Ri, N)

% Calculate regression rate in g/(cm^2*s)
% Could probably be split into its own function
Gox = (moDot*1000)./(pi*(Ri.*100).^2); % g/(cm^2*s)
rDot = (a.*Gox.^n)./1000; % m/s

% Calculate the required fuel grain length
L = (mfDot./N)./(2*pi.*Ri.*rhoF.*rDot);

end