%% Anthony Torres
% COBRA - Hybrid Dev
% calcThroatArea.m
%
% Created: 2/9/17
% Modified: 2/9/17
%
% Purpose: Calculates the required throat area for a given starting area
% and assuming M=1 at the throat
% Uses Braeunig (1.26) http://www.braeunig.us/space/propuls.htm
%
% Limitations:
%   
%
% Inputs:
%   mDot --- choked/max mass flow rate (kg/s)
%
%   Pt --- pressure at throat (Pa)
%
%   Tt --- temperature at throat (K)
%
%   Mm --- average molar mass of gas flowing through throat (kg/mol)
%
%   gamma --- specific heat ratio (unitless)
%
% Outputs:
%   At --- throat area (m^2)
%

function[At, Dt] = calcThroatArea(mDot, Pt, Tt, Mm, gamma)
% R --- universal gas constant ( 8.3145 J/(mol*K) )
R = 8.3145;

% Calculate area
At = (mDot./Pt).*sqrt((R.*Tt)./(gamma.*Mm));

% Calculate diameter
Dt = sqrt(At./pi);

end