function [T] = CO2_rho h 2T(rho, h, T guess)
2
3 % Purpose:
4 %??????????????????????????????????????????????????????????????????????????
5 % This function places a nonlinear solver around the CO2Props function in
6 % order to change the independant variable from T & rho to rho & h.
7 %==========================================================================
8 % Inputs:
9 %??????????????????????????????????????????????????????????????????????????
10 % rho ? Density, kg/m3
11 % h ? Enthalpy, kJ/kg
12 %==========================================================================
13 % Outputs:
14 %??????????????????????????????????????????????????????????????????????????
15 % T ? Temperature, K
16 %==========================================================================
17
18 rho Known = rho; % Save rho to a more intuitive name
19 h Known = h; % Save h to a more intuitive name
20 %T guess = 300; % Guess value for T for lsqnonlin to use
21
22 % Create a function for lsqnonlin to solve T(rho,h)
23 pFunc = @(T Unknown) getfield(CO2Props(T Unknown,rho Known),'h')?h Known;
24 %pFunc = @(T Unknown) ...
getfield(CO2PropsNIST(T Unknown,rho Known),'h')?h Known;
25 % Sinse T Unknown is not pre defined in pFunc, lsqnonlin will find a T for
26 % rho Known and h Known
70
27 T = lsqnonlin(pFunc,T guess,0,inf,optimset('Display','off'));
28
29 end