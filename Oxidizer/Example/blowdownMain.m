clear all; clc;
close all;
% Initial Conditions
%==========================================================================
% Tank Properties
%??????????????????????????????????????????????????????????????????????????
V = .0498; % Tank volume, m3
M_o = 12.76; % Initial tank fluid mass, kg
T1_o = 274.25; % Initial tank fluid temperature, K
%M o = 10.15;
%T1 o = 277.9;
%P o = 4.025; % Initial tank fluid pressure, MPa
%X o = 0.01; % Initial tank fluid quality,
R = 0.1889241; % Gas constant, kJ/(kg*K)
% Injector Properties
%??????????????????????????????????????????????????????????????????????????
d_inj = 0.178; % Injector diameter, in
Ac = (pi/4)*(d_inj*.0254)^2; % Injector cross sectional area; m3
Cd = 0.8; % Injector discharge coefficient
Pamb = 85.9e-3; % Atmospheric pressure, MPa. (4500ft)
Inj Switch = 1; % Isentropic = 1, Adiabatic = 2
% Time Iteration
%??????????????????????????????????????????????????????????????????????????
tstop = 30; % Stop time, s
dt = .5; % Step time, s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize state
%??????????????????????????????????????????????????????????????????????????
rho1_o = M_o/V;
Props1_o = CO2Props(T1_o,rho1_o);
%Props1_o = CO2PropsNIST(T1 o,rho1 o);
P1_o = Props1_o.P;
rhol_o = Props1_o.rho_l;
rhov_o = Props1_o.rho_v;
X1_o = (rhov_o/rho1_o)*(rhol_o-rho1_o)/(rhol_o-rhov_o);
h1_o = Props1_o.h;
H1_o = M_o*Props1_o.h;
st_o = Props1_o.state;

t = 0; % Column 1: Time (sec)
M = M_o; % Column 2: Tank fluid mass, M (kg)
rho1 = rho1_o; % Column 3: Tank fluid density, rho1 (kg/m�3)
T1 = T1_o; % Column 4: Tank temperature, T1 (K)
P1 = P1_o; % Column 5: Tank pressure, P1 (MPa)
X1 = X1_o; % Column 6: Tank quality, X1 ()
h1 = h1_o; % Column 7: Tank specific enthalpy, h1 (kJ/kg)
H1 = H1_o; % Column 8: Tank total enthalpy, H1 (J) ??? Units
mdot = 0; % Column 9: Tank mass flow rate, mdot (kg/s)
P2 = Pamb; % Column 10: Injector outlet pressure, P2 (Pa)
st = st_o; % Column 11: Fluid state, ?1=? Input, 0=Liq, 1=Sat, 2=Gas


State = [t, M, rho1, T1, P1, X1, h1, H1, mdot, P2, st];
State2 = [0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 0;
while t < tstop

i = i + 1;

% Exit the loop when the tank is out of fluid
%??????????????????????????????????????????????????????????????????????
if M < 0; break; end

% Initial Tank Fluid Properties
%======================================================================
%Note: this guess MUST yeild a quality less than 1 and >0
guess=[300 300]; %[T,rho]
%
%Set up function (need to match pressure and quality)
pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P1; ...
getfield(CO2Props(v(1),v(2)),'X')-X1];
%{
%Set up function (need to match pressure and quality)
pFunc = @(v) [getfield(CO2PropsNIST(v(1),v(2)),'P')?P1 ...
getfield(CO2PropsNIST(v(1),v(2)),'X')?X1];
%}
% Solve for [T,rho] of saturated but pure liquid at P1
% lsqnonlin tries to find a T and rho that makes pFunc = 0
v1 = ...
lsqnonlin(pFunc,guess,0,inf,optimset('Display','off','TolFun',1e-14));
T1 = v1(1); rho1 = v1(2);
Props1 = CO2Props(T1,rho1);
% Props1 = CO2PropsNIST(T1,rho1);
Pv1 = Props1.P; % Fluid Vapor Pressure, MPa
rhoL1 = Props1.rho_l; % Fluid liquid density, kg/m3
h1 = Props1.h; % Fluid specific enthalpy, kJ/kg
H1 = M*h1; % Fluid total enthalpy, kJ
s1 = Props1.s; % Fluid entropy, kJ/kg*K

% Propegate properties across the injector
%======================================================================

%Assume an isentropic (s1=s2) or adiabatic (h1=h2) injector , solve for properties at P2
if Inj_Switch == 1
%Set up function to match pressure and entropy
%
pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P2 ...
getfield(CO2Props(v(1),v(2)),'s')-s1];
%{
pFunc = @(v) [getfield(CO2PropsNIST(v(1),v(2)),'P')?P2 ...
getfield(CO2PropsNIST(v(1),v(2)),'s')?s1];
%}
elseif Inj_Switch == 2
%Set up function to match pressure and enthalpy
pFunc = @(v) [getfield(CO2Props(v(1),v(2)),'P')-P2 ...
getfield(CO2Props(v(1),v(2)),'h')-h1];
end
% Solve for T2 & rho2 downstream of the injector
%??????????????????????????????????????????????????????????????????????
guess=[300 300]; %[T,rho]
v2 = ...
lsqnonlin(pFunc,guess,0,inf,optimset('Display','off','TolFun',1e-14));
T2 = v2(1); rho2 = v2(2);
% Solve for the remaining properties downstream of the injector
%??????????????????????????????????????????????????????????????????????
Props2=CO2Props(T2,rho2); % Downstream fluid properties
% Props2=CO2PropsNIST(T2,rho2);
Pv2 = Props2.P; % Downstream fluid vapor pressure
h2 = Props2.h; % Downstream fluid enthalpy

% Calculate the Mass Flow
%======================================================================
% Non?Equalibrium Parameter
%k = sqrt(abs(P1?P2)/abs(Pv2?P1)); % Whitmores Equation
k = sqrt((P1-P2)/(Pv1-P2)); % Spencer and Stanfords Equation
% Weighting Coefficient

W = (1/(k+1));
% Incompressible fluid mass flow rate
mdot_inc = Ac*sqrt(2*rhoL1*(P1-P2)*1e6);
% Homogeneous Equilibrium mass flow rate
mdot_HEM = rho2*Ac*sqrt(2*(h1-h2));
% Weighted Non?Homogeneous Equilibrium (modified Stanford) mass flow rate
mdot = Cd*((1-W) * mdot_inc + W * mdot_HEM); % Shannon's Theory

% Update the upstream fluid properties for the next step
%======================================================================
M = M - mdot*dt; % Update tank fluid mass
Hdot = h1*mdot; % Enthalpy flow rate
H1 = H1 - Hdot*dt; % Update tank total enthalpy

% Calculate the new tank enthalpy and density
%?????????????????????????????????????????????????????????????????????
rho1 = M/V; % Update tank specific density
h1 = H1/M; % Update tank specific enthalpy

% Calculate the new tank temperature
%??????????????????????????????????????????????????????????????????????
% Create a function for lsqnonlin to solve T(rho,h)
pFunc = @(T_Unknown) getfield(CO2Props(T_Unknown,rho1),'h')-h1;
% pFunc = @(T Unknown) getfield(CO2PropsNIST(T Unknown,rho1),'h')?h1;
% Sinse T Unknown is not pre defined in pFunc, lsqnonlin will find a ...
% T for rho Known and h Known
T1 = lsqnonlin(pFunc,300,0,inf,optimset('Display','off','TolFun',1e-14));
% [T1] = CO2 rho h 2T(rho1,h1,300);

% Calulate the new tank pressure and quality
%??????????????????????????????????????????????????????????????????????
Props1 = CO2Props(T1,rho1);
% Props1 = CO2PropsNIST(T1,rho1);
P1 = Props1.P;
X1 = Props1.X;

st1 = Props1.state;

% Update the state
%??????????????????????????????????????????????????????????????????????
t = t + dt;

State = [State; t, M, rho1, T1, P1, X1, h1, H1, mdot, Pamb, st1];
State2 = [State2; mdot inc, mdot HEM];
end
%
% Open cold flow data for comparison
%==========================================================================
%  CF = open('Cold Flow Data\FillTest1 29 largeNoTop.mat');
%CF = open('Cold Flow Data\FillTest1 29 smallNoTop.mat');
%CF = open('Cold Flow Data\FillTest1 29 large.mat');
%CF = open('Cold Flow Data\FillTest1 29 small.mat');

% Find portion of cold flow data where the run valve is open
%??????????????????????????????????????????????????????????????????????????
%  iter = 0;
%  for loop = 1:length(CF.Valvedat)
%  if CF.Valvedat(loop,5) == 1
%  iter = iter + 1;
%  CF_t(iter) = CF.t(loop);
%  CF_Weight(iter) = CF.Weight(loop);
%  CF_PTank(iter) = CF.PRunTank(loop);
%  CF_Mdot(iter) = CF.mdotVentIn(loop);
%  CF_DeltaWeight(iter) = CF.DeltaWeightKg(loop);
%  CF_T_Tank(iter) = CF.Tank_Temp(loop);
%  CF_dP_Vent(iter) = CF.dpVenturi(loop);
%  CF_rhol_H(iter) = CF.rhol_H(loop);
%  end
%  end
%  CF_t = CF_t - CF_t(1);


% Calculate the venturi mass flow rate the helmholtz density
%??????????????????????????????????????????????????????????????????????????
d1=0.0254762;
d2=0.0078486;
A1=pi*d1^2/4;
A2=pi*d2^2/4;
C=0.980243097172394;
%  CF_Mdot_H=A2.*C.*CF_rhol_H.*sqrt(2.*(CF_dP_Vent)./(CF_rhol H.*(1-(A2/A1)^2)));
% 
%  % Calculate the mass flow rate by differntiating the tank mass
%  %??????????????????????????????????????????????????????????????????????????
%  CF Mdot diff = ?diff(CF Weight)./diff(CF t);
%  CF dt = CF t(2) ? CF t(1);
%  %CF Mdot wave = derivative cwt(CF Weight,'gaus1',16,CF dt,1);
%  CF Mdot wave = derivative dwt(CF Weight,'spl',4,CF dt,1);
% 
%  % Open cold flow data for comparison
%  %==========================================================================
%  Whitmore = open('Whitmores Output\FillTest1 29 largeNoTop.mat');
%  W t = Whitmore.t;
%  W M = Whitmore.M;
%  W rho1 = Whitmore.rho1;
%  W T1 = Whitmore.T1;
%  W P1 = Whitmore.P1;
%  W X1 = Whitmore.X1;
%  W H1 = Whitmore.H1;
%  W mdot = Whitmore.mdot;
%  W 2 = Whitmore.P2;

% Up?pack the state vector for ploting
%==========================================================================
t = State(:,1); % Column 1: Time (sec)
M = State(:,2); % Column 2: Tank fluid mass, M (kg)
rho1 = State(:,3); % Column 3: Tank fluid density, rho1 (kg/m�3)
T1 = State(:,4); % Column 4: Tank temperature, T1 (K)

P1 = State(:,5); % Column 5: Tank pressure, P1 (Pa)
X1 = State(:,6); % Column 6: Tank quality, X1 ()
h1 = State(:,7); % Column 7: Tank specific enthalpy, h1 (kJ/kg)
H1 = State(:,8); % Column 8: Tank total enthalpy, H1 (J) ??? Units
mdot = State(:,9); % Column 9: Tank mass flow rate, mdot (kg/s)
P2 = State(:,10); % Column 10: Injector outlet pressure, P2 (Pa)
st = State(:,11); % Column 11: Fluid state, ?1=? Input, 0=Liq, 1=Sat, 2=Gas

mdot_inc = State2(:,1); % Mdot Incompressible
mdot_HEM = State2(:,2); % Mdot HEM

%% Plot
%==========================================================================
% Column 2 ? Tank fluid mass
%??????????????????????????????????????????????????????????????????????????
figure(2)
plot(t,M); hold on; grid on;
%plot(W t,W M,'??g');
%  plot(CF t,CF Weight,'??r');
%plot(CF t,CF DeltaWeight,'g');
%title('Tank Mass');
xlim([0 tstop]);
xlabel('Time, s'); ylabel('Tank Fluid Mass, kg');
%legend('Model Predicted','Whitmore','Cold Flow');
legend('Model Predicted','Cold Flow');
%  saveas(2,'Tank Mass.png')

% Column 3 ? Tank fluid specific density
%??????????????????????????????????????????????????????????????????????????
figure(3);
plot(t,rho1); hold on; grid on;
plot(W_t,W_rho1,'??g');
%  plot(CF_t,CF Weight/V,'r');
title('rho1');

xlabel('Time, s'); ylabel('Density, kg/m3');
legend('Model Predicted','Whitmore','Cold Flow');
%  saveas(3,'Tank Density.png')

% Column 4 ? Tank temperature
%??????????????????????????????????????????????????????????????????????????
figure(4);
plot(t,T1); hold on; grid on;
%plot(W t,W T1,'??g');
%  plot(CF_t,CF T Tank,'??r');
%title('Temperature');
xlim([0 tstop]);
xlabel('Time, s'); ylabel('Temperature, K');
%legend('Model Predicted','Whitmore','Cold Flow');
legend('Model Predicted','Cold Flow');
%  saveas(4,'Tank Temperature.png')

% Column 5 ? Tank pressure
%??????????????????????????????????????????????????????????????????????????
figure (5)
plot(t,P1); hold on; grid on;
%plot(W t,W P1/1e3,'??g');
%  plot(CF t,CF PTank/1e6,'??r');
%title('P1');
xlim([0 tstop]);
xlabel('Time, s'); ylabel('Pressure, MPa');
legend('Model Predicted','Cold Flow');
%legend('Model Predicted','Whitmore','Cold Flow');
%  saveas(5,'Tank Pressure.png')

% Column 6 ? Tank quality
%??????????????????????????????????????????????????????????????????????????
figure (6)
plot(t,X1); hold on; grid on;
plot(W_t,W_X1,'??g');

ylim([0 1]);
title('X1');
xlabel('Time, s'); ylabel('Quality');
legend('Model Predicted','Whitmore');
%  saveas(6,'Tank Quality.png')

% Coumn 7 ? Tank specific enthalpy
%??????????????????????????????????????????????????????????????????????????
figure (7)
plot(t,h1); hold on; grid on;
title('h1');
xlabel('Time, s'); ylabel('Specific Enthalpy, kJ/kg?K');
legend('Model Predicted');
%  saveas(7,'Tank Specific Enthalpy.png')

% Coumn 8 ? Tank total enthalpy
%??????????????????????????????????????????????????????????????????????????
figure (8)
plot(t,H1); hold on; grid on;
plot(W_t,W_H1,'??g');
title('H1'); xlabel('time, s'); ylabel('Enthalpy');
legend('Model Predicted','Whitmore');
%  saveas(8,'Tank Total Enthalpy.png')

% Column 9 ? Mass flow rate
%??????????????????????????????????????????????????????????????????????????
figure(9);
plot(t,mdot,'b'); hold on; grid on;
%  plot(CF t,CF Mdot,'??r')
%  plot(CF t,CF Mdot H,'?.b')
xlabel('Time, s'); ylabel('Mass Flow Rate, kg/s');
xlim([0 tstop]); %ylim([0 3]);
legend('Model Predicted','Venturi?Curve Fit Density','Venturi?Helmholtz Density');
%  saveas(9,'Tank Mass Flow Rate Venturi.png')


figure(90);
plot(t,mdot,'b'); hold on; grid on;
%plot(W t,W mdot,'??g');
%  plot(CF t,CF Mdot,'??r')
%  plot(CF t,?CF Mdot wave,'?.b')
plot(t,mdot_inc,':b');
plot(t,mdot_HEM,'??bo');
xlabel('Time, s'); ylabel('Mass Flow Rate, kg/s');
xlim([0 tstop]); ylim([0 1.6]);
legend('Model Predicted','Cold Flow?Venturi','Cold Flow?Differentiated Load','Incompressible Model','HEM Model');
%  saveas(90,'Tank Mass Flow Rate.png')

%  figure(900);
%  plot(CF t(2:end),CF Mdot diff); hold on; grid on;
%  plot(CF t,?CF Mdot wave,'??r')
%  xlabel('Time, s'); ylabel('Mass Flow Rate, kg/s');
%  xlim([0 tstop]); ylim([?2 4]);
%  legend('Cold Flow?dM/dt?diff','Cold Flow?dM/dt?wavelet');
%  saveas(900,'Tank Mass Flow Rate Diff.png')

% Column 11 ? Tank fluid state
%??????????????????????????????????????????????????????????????????????????
figure(11);
plot(t,st); grid on;
title('Fluid State');
xlabel('Time, s'); ylabel('Fluid State');
ylim([2 3]);
legend('Model Predicted');
%  saveas(11,'Tank Fluid State.png')