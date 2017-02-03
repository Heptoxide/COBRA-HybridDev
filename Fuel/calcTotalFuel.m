%% Anthony Torres
% COBRA - Hybrid Dev
% calcTotalFuel.m
%
% Created: 2/3/17
% Modified: 2/3/17
%
% Purpose: Calculates the required total fuel mass using (15-17) of RPE
%
% Limitations:
%   Only works for circular port geometries
%
% Inputs:
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
%   L --- length of fuel req, see calcFuelLength.m (m)
%
%   t --- burn time (s)
%
% Outputs:
%   mF --- necessary mass of fuel
%

function[mF] = calcTotalFuel(moDot, rhoF, a, n, Ri, N, L, t)

% Calculate the necessary fuel mass given the input parameters
% NOTE: a is divided by 1000 to convert from mm/s to m/s
mF = pi*N*rhoF*L.*(((a./1000).*(2*n+1).*(moDot./(pi.*N)).^n.*t + Ri.^(2*n+1)).^(2./(2*n+1)) - Ri.^2);

end