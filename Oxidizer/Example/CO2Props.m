function [Props] = CO2Props(T,rho)
%% INTRODUCTION
%==========================================================================
% Purpose:
% Returns a sturcture containing the thermodynamic properties for carbon
% dioxide (CO2) for a given density (kg/m^3) and temperature (K). Valid
% temperature and pressure ranges are:
% 216 K <= T <= 1100K       0 MPa <= P <= MPa
%
% The properties in the structure are:
%
% P X s u cv cp h c
% rho l s l u l cv l cp l h l c l
% rho v s v u v cv v cp v h v c v
%
% Inputs:
% T - Temperature, K
% rho - Density, kg/m3
% Outputs:
% P - Pressure, MPa
% X - Quality, (Vapor Mass/Total Fluid Mass)
% s - Specific Entropy, kJ/(kg*K)
% u - Specific Internal Energy, kJ/kg
% cv - Specific Heat at Constant Volume, kJ/(kg*K)
% cp - Specific Heat at Constant Pressure, kJ/(kg*K)
% h - Specific Enthalpy, kJ/kg
% c - Speed of Sound, m/s
% rho - Density, kg/m3
% state - -1 = Negative Input, 0 = Liquid, 1 = Saturated, 2 = Gas
% l - Liquid designator
% v - Vapor designator
%--------------------------------------------------------------------------
% Revision History:

% Written for a CO2 Blowdown model developed at Utah State University by
% Matthew Wilson
% 4130 Old Main Hill
% Logan, UT 84322-4130
% Recommented and checked by:
% Brian Solomon
%--------------------------------------------------------------------------
% Based upon the Helmholtz Energy based equations of state described by
% Span, R. and Wagner, W. in
% "A New Equation of State for Carbon Dioxide Covering the Fluid Region
% from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa"
% Journal of Physical and Chemical Reference Data
% Vol 25, No. 6, 1996. Pp 1509-1596
%--------------------------------------------------------------------------

%% CALCULATE THE SATURATION PROPERTIES
%==========================================================================
% CO2 Constants
%--------------------------------------------------------------------------
R = 0.1889241; % Gas constant, kJ/(kg*K)
Tt = 216.592; % Triple point temperature, K (Eq. 3.1)
Pt = 0.51795; % Triple point pressure, MPa (Eq. 3.2)
Tc = 304.1282; % Critical temperature, K (Eq. 3.3)
Pc = 7.3773; % Critical pressure, MPa (Eq. 3.4)
rhoc = 467.6; % Critical density, kg/m3 (Eq. 3.5)
% d = rho/rhoc; % Reduced density,  = phi/phi critical, phi - mass density
t = Tc/T; % Inverse reduced temperature, tau = T critical/T, T - temperature

% Coefficients
%--------------------------------------------------------------------------
% Melting pressure coefficients, Section 3.3
a_m = [1955.5390 2055.4593]';
% Sublimation pressure coefficients, Section 3.4

a_s = [-14.740846 2.4327015 -5.3061778]';
% Vapor pressure coefficients, Section 3.5
a_p = [-7.0602087 1.9391218 -1.6463597 -3.2995634]';
t_p = [1.0 1.5 2.0 4.0]';
% Saturated liquid density coeffiecients, Section 3.6
a_l = [1.9245108 -0.62385555 -0.32731127 0.39245142]';
t_l = [0.34 0.5 10/6 11/6]';
% Saturated vapor density coeffiecients, Section 3.7
a_v = [-1.7074879 -0.82274670 -4.6008549 -10.111178 -29.742252]';
t_v = [0.34 0.5 1.0 7/3 14/3]';

% Phase Property Calculations
%--------------------------------------------------------------------------
% P melt = Pt* (1 + a m(1)*(T/Tt-1) + a m(2)*(T/Tt-1)^2); % Melting ...
% Pressure, Eq. 3.10
% P sub = Pt*exp(Tt/T*(a s(1)*(1-T/Tt) + ... % Sublimation ...
% Presssure, Eq. 3.12
% a s(2)*(1-T/Tt)^1.9 + ...
% a s(3)*(1-T/Tt)^2.9));
P_sat = Pc * exp(Tc/T*sum(a_p.*(1-T/Tc).^t_p)); % Vapor Pressure, Eq. 3.13
rho_l = rhoc * exp(sum(a_l.*(1-T/Tc).^t_l)); % Saturated Liquid Density, Eq. 3.14
rho_v = rhoc * exp(sum(a_v.*(1-T/Tc).^t_v)); % Saturated Vapor Density, Eq. 3.15

%% CALCULATE THE PROPERTIES AT THE GIVEN TEMPERATURE AND DENSITY
%==========================================================================
% Account for negative density or temperature
%--------------------------------------------------------------------------
if rho < 0 || T < 0
    X = NaN;
    P = NaN;
    s_v = NaN; u_v = NaN; cp_v = NaN;
    cv_v = NaN; h_v = NaN; c_v = NaN;
    s_l = NaN; u_l = NaN; cp_l = NaN;
    cv_l = NaN; h_l = NaN; c_l = NaN;
    
    s = s_v*X + s_l*(1-X);
    u = u_v*X + u_l*(1-X);
    h = h_v*X + h_l*(1-X);
    cp = cp_v*X + cp_l*(1-X);
    cv = cv_v*X + cv_l*(1-X);
    c = NaN;
    state = -1;
% GAS
%--------------------------------------------------------------------------
elseif rho < rho_v || imag(rho_l)
    [P, s, u, cp, cv, h, c] = Helmholtz(t,rho/rhoc);
    X = 1;
    s_v = s; u_v = u; cp_v = cp;
    cv_v = cv; h_v = h; c_v = c;
    s_l = NaN; u_l = NaN; cp_l = NaN;
    cv_l = NaN; h_l = NaN; c_l = NaN;
    state = 2;
    %P = P/1E6;
% LIQUID
%--------------------------------------------------------------------------
elseif rho > rho_l
    [P, s, u, cp, cv, h, c] = Helmholtz(t,rho/rhoc);
    X = 0;
    s_l = s; u_l = u; cp_l = cp;
    cv_l = cv; h_l = h; c_l = c;
    s_v = NaN; u_v = NaN; cp_v = NaN;
    cv_v = NaN; h_v = NaN; c_v = NaN;
    state = 0;
    %P = P/1E6;
% MELTING---
%--------------------------------------------------------------------------
%
% SATURATED
%--------------------------------------------------------------------------

else
    X = rho_v*(rho_l - rho)/(rho*(rho_l-rho_v));
    [P, s_l, u_l, cp_l, cv_l, h_l, c_l] = Helmholtz(t,rho_l/rhoc);
    [P, s_v, u_v, cp_v, cv_v, h_v, c_v] = Helmholtz(t,rho_v/rhoc);
    s = s_v*X + s_l*(1-X);
    u = u_v*X + u_l*(1-X);
    h = h_v*X + h_l*(1-X);
    cp = cp_v*X + cp_l*(1-X);
    cv = cv_v*X + cv_l*(1-X);

    c = NaN;
    P = P_sat;
    state = 1;
end

%% CREATE THE OUTPUT STRUCTURE
%==========================================================================
Props.P = P; % Pressure
Props.X = X; % Quality
Props.s = s;            Props.s_l = s_l;        Props.s_v = s_v;    % Entropy
Props.u = u;            Props.u_l = u_l;        Props.u_v = u_v;    % Internal Energy
Props.cv = cv;          Props.cv_l = cv_l;      Props.cv_v = cv_v;  % Specific Heat at Constant Volume
Props.cp = cp;          Props.cp_l = cp_l;      Props.cp_v = cp_v;  % Specific Heat at Constant Pressure
Props.h = h;            Props.h_l = h_l;        Props.h_v = h_v;    % Enthalpy
Props.c = c;            Props.c_l = c_l;        Props.c_v = c_v;    % Speed of Sound
Props.rho_l = rho_l;    Props.rho_v = rho_v;                        % Density
Props.state = state;                                                % State
%Props.gma l = cp l/cv l;
%Props.gma v = cp v/cv v;
end

%##########################################################################
function [P, s, u, cp, cv, h, c] = Helmholtz(t,d)
%
%% INTORDUCTION
%==========================================================================
% Calculates CO2 properties using the residual and ideal gas portions of
% the Helmholtz Energy at a given temperature and density normalized to the
% critical point temperature and density.
%--------------------------------------------------------------------------
% Based upon the Helmholtz Energy based equations of state described by
% Span, R. and Wagner, W. in
% "A New Equation of State for Carbon Dioxide Covering the Fluid Region
% from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa"
% Journal of Physical and Chemical Reference Data
% Vol 25, No. 6, 1996. Pp 1509-1596
%--------------------------------------------------------------------------
% Input:
% t - Inverse Reduced Temperature, tau = T critical/T, T - Temperature
% d - Reduced Density,  = phi/phi critical, phi - Mass Density
%--------------------------------------------------------------------------
% Output:
% P - Pressure, Pa
% s - Entropy, kJ/(kg*K)
% u - Internal energy, kJ/kg
% cp - Isochoric heat capacity, kJ/(kg*K)
% cv - Isobaric heat capacity, kJ/(kg*K)
% h - Entropy, kJ/kg
% c - Speed of sound, m/s
%==========================================================================

%% RESIDUAL PART OF HELMHOLTZ
%==========================================================================
% CONSTANT PARAMETER DEFINITION -
% Page 1544, Table 31. Coefficients and exponents of Eq. (6.5)
%--------------------------------------------------------------------------
n7 = [ 0.38856823203161E0
        0.29385475942740E1
        -0.55867188534934E1
        -0.76753199592477E0
        0.31729005580416E0
        0.54803315897767E0
        0.12279411220335E0];
d7 = [1 1 1 1 2 2 3]';
t7 = [0.00 0.75 1.00 2.00 0.75 2.00 0.75]';
%--------------------------------------------------------------------------
n34 = [ 0.21658961543220E1
        0.15841735109724E1
       -0.23132705405503E0
        0.58116916431436E-1
%         
       -0.55369137205382E0
        0.48946615909422E0
       -0.24275739843501E-1
        0.62494790501678E-1
%         
       -0.12175860225246E0
       -0.37055685270086E0
       -0.16775879700426E-1
       -0.11960736637987E0
%         
       -0.45619362508778E-1
        0.35612789270346E-1
       -0.74427727132052E-2
       -0.17395704902432E-2
%         
       -0.21810121289527E-1
        0.24332166559236E-1
       -0.37440133423463E-1
        0.14338715756878E0
%         
       -0.13491969083286E0
       -0.23151225053480E-1
        0.12363125492901E-1
        0.21058321972940E-2
%         
       -0.33958519026368E-3
        0.55993651771592E-2
       -0.30335118055646E-3];

d34 = [1 2 4 5 5 5 6 6 6 1 1 4 4 4 7 8 2 3 3 5 5 6 7 8 10 4 8]';
t34 = [1.5000 1.5000 2.5000 0.0000 1.5000 2.0000 0.0000 1.0000...
    2.0000 3.0000 6.0000 3.0000 6.0000 8.0000 6.0000 0.0000 7.0000 12.0000 ...
    16.0000 22.0000 24.0000 16.0000 24.0000 8.0000 2.0000 28.0000 14.0000]';
c34 = [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 4 4 4 4 4 4 5 6]';
%--------------------------------------------------------------------------
n39 = [-0.21365488688320E3
        0.26641569149272E5
       -0.24027212204557E5
       -0.28341603423999E3
        0.21247284400179E3];
d39 = [2 2 2 3 3]';
t39 = [1.00 0.00 1.00 3.00 3.00]';
a39 = [25 25 25 15 20]';
bt39 = [325 300 300 275 275]';
y39 = [1.16 1.19 1.19 1.25 1.22]';
e39 = [1.00 1.00 1.00 1.00 1.00]';
%--------------------------------------------------------------------------
n42 = [-0.66642276540751E0
    0.72608632349897E0
    0.55068668612842E-1];
a42 = [3.500 3.500 3.000]';
b42 = [0.875 0.925 0.875]';
bt42 = [0.300 0.300 0.300]';
A42 = [0.700 0.700 0.700]';
B42 = [0.3 0.3 1.0]';
C42 = [10.0 10.0 12.5]';
D42 = [275 275 275]';

R = 0.1889241; % Gas Constant, kJ/(kg*K)
Tc = 304.1282; % Critical Temperature, K
rhoc = 467.6; % Critical Density, kg/m3
%--------------------------------------------------------------------------

% NONDIMENSIONALIZED HELMHOLTZ FREE ENERGY

% RESIDUAL PORTION OF THE HELMHOLTZ ENERGY
%--------------------------------------------------------------------------
tta = (1-t) + A42 .* ((d-1).^2).^(0.5./bt42);
Del = tta.^2 + B42.*((d-1).^2).^a42;
Psi = exp(-C42.*(d-1).^2 - D42.*(t-1).^2);

Psi_d = -2 .* C42 .* (d-1) .* Psi;
Psi_t = -2 .* D42 .* (t-1) .* Psi;
Del_d = (d-1).*(A42 .* tta .* 2./bt42 .* ...
    ((d-1).^2).^(.5./bt42-1) + ...
    2.* B42 .* a42 .* ((d-1).^2).^(a42-1));

Psi_dd = 2.*C42.*Psi.*(2.*C42.*(d-1).^2 - 1);
Psi_tt = 2.*D42.*Psi.*(2.*D42.*(t-1).^2 - 1);
Del_dd = 1./(d-1).*Del_d + (d-1).^2 ...
    .* (4.*B42.*a42.*(a42-1).*((d-1).^2).^(a42-2) + ...
    2.*A42.^2.*(1./bt42).^2.*(((d-1).^2).^(1./(2.*bt42)-1)).^2 + ...
    A42.*tta.*4./bt42.*(.5./bt42-1).*((d-1).^2).^(.5./bt42-2));

Psi_dt = 4.*C42.*D42.*(d-1).*(t-1).*Psi;

% Residual part of the Helmholtz energy, Eq. 6.5
%--------------------------------------------------------------------------
phir = sum(n7 .* d.^d7 .* t.^t7) ...
    + sum(n34 .* d.^d34 .* t.^t34 .* exp(-d.^c34)) ...
    + sum(n39 .* d.^d39 .* t.^t39 .* ...
        exp(-a39.*(d - e39).^2 - bt39.*(t - y39).^2)) ...
    + sum(n42 .* Del.^b42 .* d .* Psi);

phir_d = sum(n7 .* d7 .* d.^(d7-1) .* t.^t7) ...
    + sum(n34 .* d.^(d34-1) .* t.^t34 .* ...
        exp(-d.^c34) .* (d34 - c34.* d.^c34)) ...
    + sum(n39 .* d.^d39 .* t.^t39 .* ...
        exp(-a39.*(d - e39).^2 - bt39.*(t - y39).^2) .* ...
        (d39./d - 2.*a39.*(d-e39))) ...
    + sum(n42 .* (Del.^b42 .* (Psi + d.*Psi_d) + ...
        b42 .* Del.^(b42-1) .* Del_d .* d .* Psi));

phir_dd = sum(n7 .* d7 .* (d7 - 1) .* d.^(d7-2) .* t.^t7) ...
    + sum(n34 .* exp(-d.^c34) .* d.^(d34-2) .* t.^t34 .* ...
        ((d34 - c34.* d.^c34).*(d34 -1-c34.* d.^c34) - c34.^2 .* ...
            d.^c34)) ...
    + sum(n39 .* t.^t39 .* ...
        exp(-a39.*(d - e39).^2 - bt39.*(t - y39).^2) .* ...
        (-2.*a39.*d.^d39 + 4.*a39.^2.*d.^d39.*(d-e39).^2 - ...
        4.*d39.*a39.*d.^(d39-1).*(d-e39) + d39.*(d39-1).*d.^(d39-2))) ...
    + sum(n42 .* (Del.^b42 .* (2.*Psi_d + d*Psi_dd) + ...
        2.*b42 .* Del.^(b42-1) .* Del_d .*(Psi + d.*Psi_d) + ...
        b42.*d.*Psi.*...
            (Del.^(b42-1).*Del_dd + (b42-1).*Del.^(b42-2).*Del_d.^2)));

phir_t = sum(n7 .* t7 .* d.^d7 .* t.^(t7-1)) ...
    + sum(n34 .* d.^d34 .* t34 .* t.^(t34-1) .* exp(-d.^c34)) ...
    + sum(n39 .* d.^d39 .* t.^t39 .* ...
        exp(-a39.*(d - e39).^2 - bt39.*(t - y39).^2) .* ...
        (t39./t - 2.*bt39.*(t-y39))) ...
    + sum(n42 .*d .* (-2.*tta.*b42.*Del.^(b42-1).*Psi + ...
        Del.^b42.*Psi_t));

phir_tt = sum(n7 .* t7 .* (t7-1) .* d.^d7 .* t.^(t7-2)) ...
    + sum(n34 .* t34 .* (t34-1) .* d.^d34 .* t.^(t34-2) .* ...
        exp(-d.^c34)) ...
    + sum(n39 .* d.^d39 .* t.^t39 .* ...
        exp(-a39.*(d - e39).^2 - bt39.*(t - y39).^2) .* ...
        ((t39./t - 2.*bt39.*(t-y39)).^2 - t39./(t.^2)-2.*bt39)) ...
    + sum(n42 .*d .* (Psi.*(2.*b42.*Del.^(b42-1) + ...
        4.*tta.^2.*b42.*(b42-1).*Del.^(b42-2)) - ...
        4.*tta.*b42.*Del.^(b42-1).*Psi_t + Del.^b42.*Psi_tt));

phir_dt = sum(n7 .* t7 .* d7 .* d.^(d7-1) .* t.^(t7-1)) ...
    + sum(n34 .* exp(-d.^c34) .* t34 .* d.^(d34-1) .* t.^(t34-1) .* ...
        (d34 - c34.* d.^c34)) ...
    + sum(n39 .* d.^d39 .* t.^t39 .* ...
        exp(-a39.*(d - e39).^2 - bt39.*(t - y39).^2) .* ...
        (t39./t - 2.*bt39.*(t-y39)).*(d39./d - 2.*a39.*(d-e39))) ...
    + sum(n42.* ( Del.^b42.*(Psi_t + d.*Psi_dt) + ...
        d.*b42.*Del.^(b42-1).*Del_d.*Psi_t - ...
        2.*tta.*b42.*Del.^(b42-1).*(Psi + d.*Psi_d) + ...
        d.*Psi.*(-A42.*b42.*0.5./bt42.*Del.^(b42-1).*(d-1).*((d-1).^2).^(0.5./bt42-1) ...
    - ...
        2.*tta.*b42.*(b42-1).*Del.^(b42-2).*Del_d)));

%% IDEAL GAS HELMHOLTZ
%==========================================================================
% CONSTANT PARAMETER DEFINITION -
% Page 1540, Table 27. Coefficients of the correlation equations, Eq. (6.2)
% and Eq. (6.3)
%--------------------------------------------------------------------------
ao3 = [8.37304456 -3.70454304 2.50000000]';
ao8 = [1.99427042
    0.62105248
    0.41195293
    1.04028922
    0.08327678];
ttao8 = [3.15163
    6.11190
    6.77708
    11.32384
    27.08792];

% Ideal Gas part of the Helmholtz energy, Eq. 6.3
%--------------------------------------------------------------------------
phio = log(d) + ao3(1) + ao3(2)*t + ao3(3)*log(t)...
    + sum(ao8.*log( 1-exp(-ttao8.*t) ));
% phio d = 1/d;

% phio dd = -1/d^2;
% phio dt = 0;

phio_t = ao3(2) + ao3(3)/t + sum(ao8.*ttao8.*( (1-exp(-ttao8.*t)).^(-1) ...
-1 ));

phio_tt = -ao3(3)/(t^2) - sum(ao8 .* ttao8.^2 .* exp(-ttao8.*t) .* ...
( 1-exp(-ttao8.*t) ).^(-2) );

%% THERMODYNAMIC PROPERTY CALCULATIONS
%==========================================================================
% Page 1517, Table 3. Relations of thermodynamic properties to the
% dimensionless Helmholtz function phi consisting of phi ideal gas and
% phi residual
%--------------------------------------------------------------------------
% Pressure
P = 1000*d*rhoc*R*Tc/t*(1+d*phir_d);
% Entropy
s = R*(t*(phio_t + phir_t) - phio - phir);
% Internal energy
u = R*Tc/t*t*(phio_t+phir_t);
% Isochoric heat capacity
cv = -R*t^2*(phio_tt+phir_tt);
% Isobaric heat capacity
cp = R*(-t^2*(phio_tt+phir_tt) + ...
    (1+d*phir_d-d*t*phir_dt)^2/(1+2*d*phir_d+d^2*phir_dd));
% Enthalpy
h = R*Tc/t*(1+t*(phio_t+phir_t)+d*phir_d);
% Speed of sound
c = sqrt( R*Tc/t*1000*( 1 + 2*d*phir_d + d^2*phir_dd - ...
    ((1 + d*phir_d - d*t*phir_dt)^2/(t^2*(phio_tt+phir_tt)))) );
% Joule-Thompson coefficient
mu = -1/(R*d*rhoc)*(d*phir_d + d^2*phir_dd + d*t*phir_dt) ...
    / ((1+d*phir_d-d*t*phir_dt)^2 - ...
    t^2*(phio_tt+phir_tt)*(1+2*d*phir_d+d^2*phir_dd));

% Fugacity
%psi = exp(phir+d*phir d-log(1+d*phir d));
% 2nd Virial
%B = lim(phir d)/rhoc as d->0
% 3rd Virial
%C = lim(phir dd)/rhoc^2 as d->0

%% Convert the energy terms from J/gram to J/kg
%=========================================================================
%u = u*1000; s = s*1000; h = h*1000;
%
%% Reference the energy terms to the same reference as the NIST tables
%=========================================================================
u_NIST = 506.8006990081582;
s_NIST = 2.739101354918014;
h_NIST = 506.77604151182834;

% Shift the energy terms to match the NIST reference point
%--------------------------------------------------------------------------
u = u + u_NIST;
s = s + s_NIST;
h = h + h_NIST;
%
end