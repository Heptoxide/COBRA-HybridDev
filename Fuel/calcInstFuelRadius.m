%% Anthony Torres
% COBRA - Hybrid Dev
% calcInstFuelRadius.m
%
% Created: 2/3/17
% Modified: 2/3/17
%
% Purpose: Changed to calculate based on (15-14) of RPE
%
%   Calculates the required outer radius of a fuel grain with SINGLE
%   circular port geometry based on given fuel density, length, and mass
%   constrains
%
% Limitations:
%   Only works for SINGLE circular port geometry
%
% To Do:
%   Expand to multiple circular geometry
%
%
% Inputs:
%
%   moDot --- mass flow rate of oxidizer (kg/s)
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
%   t --- burn time (s)
%
% Outputs:
%   outerRadius --- outer radius of fuel grain (m)
%

function[radius] = calcInstFuelRadius(moDot, a, n, Ri, N, t)

% Calculate outer radius based on (15-14) of RPE
% NOTE: a is divided by 1000 to convert from mm/s to m/s
radius = ((a./1000).*(2*n+1).*(moDot./(pi.*N)).^n .* t + Ri.^(2.*n+1)).^(1/(2.*n+1));


end