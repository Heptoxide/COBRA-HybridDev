clc; clear all;

 % Create a temperature and density vector to calculate CO2 properties at
 %??????????????????????????????????????????????????????????????????????????
 T = linspace(216.59,300,50); % Temperature vector, K (?27.67 to 116.33F)
 %T = [T linspace(300,304.1,15)];
 T = [T linspace(300,307,50)];
 %rho = linspace(200,1100); % Density vector, kg/m3
 rho = 500;

 % Load the CO2 properties created by the NIST Webbook
 %??????????????????????????????????????????????????????????????????????????
 [labels,x,y] = readColData('SatCO2Props.txt',25);
 n = 3;
 T NIST = downsample(x,n);
 P NIST = downsample(y(:,1),n);
 rho l NIST = downsample(y(:,2),n);
 rho v NIST = downsample(y(:,14),n);

 u l NIST = downsample(y(:,4),n);
 u v NIST = downsample(y(:,16),n);
 h l NIST = downsample(y(:,5),n);
 h v NIST = downsample(y(:,17),n);
 s l NIST = downsample(y(:,6),n);
 s v NIST = downsample(y(:,18),n);
 cv l NIST = downsample(y(:,7),n);
 cv v NIST = downsample(y(:,19),n);
 cp l NIST = downsample(y(:,8),n);
 cp v NIST = downsample(y(:,20),n);
 c l NIST = downsample(y(:,9),n);
 c v NIST = downsample(y(:,21),n);

 % Loop through temperature and density vectors to populate the CO2 properties
 %??????????????????????????????????????????????????????????????????????????
 for i = 1:length(T)
 for j = 1:length(rho)
 Props = CO2Props(T(i),rho(j));
 P(i,j) = Props.P; % Pressure
 X(i,j) = Props.X; % Quality
 s(i,j) = Props.s; % Entropy
 s l(i,j) = Props.s l;
 s v(i,j) = Props.s v;
 u(i,j) = Props.u; % Internal Energy
 u l(i,j) = Props.u l;
 u v(i,j) = Props.u v;
 cv(i,j) = Props.cv; % Specific Heat at Constant Volume
 cv l(i,j) = Props.cv l;
 cv v(i,j) = Props.cv v;
 cp(i,j) = Props.cp; % Specific Heat at Constant Pressure
 cp l(i,j) = Props.cp l;
 cp v(i,j) = Props.cp v;
 h(i,j) = Props.h; % Enthalpy
 h l(i,j) = Props.h l;
 h v(i,j) = Props.h v;

 c(i,j) = Props.c; % Speed of Sound
 c l(i,j) = Props.c l;
 c v(i,j) = Props.c v;
 rho l(i,j) = Props.rho l; % Density
 rho v(i,j) = Props.rho v;
 state(i,j) = Props.state; % State
 end
 end

 %Shift energy terms so that the reference matches
 %??????????????????????????????????????????????????????????????????????????
 u diff = u l NIST(1) ? u l(1)
 s diff = s l NIST(1) ? s l(1)
 h diff = h l NIST(1) ? h l(1)

 % Change the independent variables to temperature and pressure.
 %??????????????????????????????????????????????????????????????????????????
 %[rhoArr gmaArr] = CO2tp(T NIST, P NIST);

 % Parse mosster data CO2 density
 %??????????????????????????????????????????????????????????????????????????
 rho parse=?0.136928567045648.*T.�2+68.960274348176327.*T?7677.810285569415;

 % Plot the pressure vs temperature
 %??????????????????????????????????????????????????????????????????????????
 figure(1);
 plot(T,P(:,1),'b'); hold on; grid on;
 plot(T NIST,P NIST,'ob');
 %title('CO2 Pressure');
 legend('Helmholtz','NIST','Location','SE');
 xlabel('Temperature, K'); ylabel('Pressure, MPa');
 saveas(1,'CO2 Pressure.png')
 %saveas(1,'LaTeX Plots/CO2 Pressure.eps','epsc')

 % Plot the density vs temperature

 %??????????????????????????????????????????????????????????????????????????
 figure(2);
 plot(T,rho l(:,1),'b'); hold on; grid on;
 plot(T NIST,rho l NIST,'ob');
 plot(T,rho v(:,1),'??r')
 plot(T NIST,rho v NIST,'sr')
 %plot(T,rho parse);
 %title('CO2 Density');
 legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST');
 %legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST','Monster');
 xlabel('Temperature, K'); ylabel('Density, kg/m3');
 saveas(2,'CO2 Density.png')
 %saveas(2,'LaTeX Plots/CO2 Density.eps','epsc')

 % Plot the internal energy vs temperature
 %??????????????????????????????????????????????????????????????????????????
 figure(3);
 plot(T,u l(:,1),'b'); hold on; grid on;
 plot(T NIST,u l NIST,'ob');
 plot(T,u v(:,1),'??r')
 plot(T NIST,u v NIST,'sr')
 title('CO2 Internal Energy');
 legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST');
 xlabel('Temperature, K'); ylabel('Internal Energy, kJ/kg');
 saveas(3,'CO2 Internal Energy.png')
 %saveas(3,'LaTeX Plots/CO2 Internal Energy.eps','epsc')

 % Plot the enthalpy vs temperature
 %??????????????????????????????????????????????????????????????????????????
 figure(4);
 plot(T,h l(:,1),'b');hold on; grid on;
 plot(T NIST,h l NIST,'ob');
 plot(T,h v(:,1),'??r')
 plot(T NIST,h v NIST,'sr')
 %title('CO2 Enthalpy');

 legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST','Location','SE');
 %legend('Liquid?NIST','Vapor?NIST');
 xlabel('Temperature, K'); ylabel('Enthalpy, kJ/kg');
 saveas(4,'CO2 Enthalpy.png')
 %saveas(4,'LaTeX Plots/CO2 Enthalpy.eps','epsc')

 % Plot the entropy vs temperature
 %??????????????????????????????????????????????????????????????????????????
 figure(5);
 plot(T,s l(:,1),'b'); hold on; grid on;
 plot(T NIST,s l NIST,'ob');
 plot(T,s v(:,1),'??r')
 plot(T NIST,s v NIST,'sr')
 %title('CO2 Entropy');
 legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST');
 xlabel('Temperature, K'); ylabel('Entropy, J/g*K');
 saveas(5,'CO2 Entropy.png')
 %saveas(5,'LaTex Plots/CO2 Entropy.eps','epsc')

 % Plot the specific heat at constant volume vs temperature
 %??????????????????????????????????????????????????????????????????????????
 figure(6);
 plot(T,cv l(:,1),'b'); hold on; grid on;
 plot(T NIST,cv l NIST,'ob');
 plot(T,cv v(:,1),'??r')
 plot(T NIST,cv v NIST,'sr')
 title('CO2 Specific Heat at Constant Volume');
 legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST');
 xlabel('Temperature, K'); ylabel('Specific Heat at Constant Volume, J/g*K');
 saveas(6,'CO2 Specific Heat at Constant Volume.png')
 %saveas(6,'LaTeX Plots/CO2 Specific Heat at Constant Volume.eps','epsc')

 % Plot the specific heat at constant pressure vs temperature
 %??????????????????????????????????????????????????????????????????????????
 figure(7);

 plot(T,cp l(:,1),'b'); hold on; grid on;
 plot(T NIST,cp l NIST,'ob');
 plot(T,cp v(:,1),'??r')
 plot(T NIST,cp v NIST,'sr')
 title('CO2 Specific Heat at Constant Pressure');
 legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST');
 xlabel('Temperature, K'); ylabel('Specific Heat at Constant Pressure, ...
J/g*K');
 saveas(7,'CO2 Specific Heat at Constant Pressure.png')
 %saveas(7,'LaTeX Plots/CO2 Specific Heat at Constant Pressure.eps','epsc')


 % Plot the speed of sound vs temperature
 %??????????????????????????????????????????????????????????????????????????
 figure(8);
 plot(T,c l(:,1),'b'); hold on; grid on;
 plot(T NIST,c l NIST,'ob');
 plot(T,c v(:,1),'??r')
 plot(T NIST,c v NIST,'sr')
 title('CO2 Speed of Sound');
 legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST');
 xlabel('Temperature, K'); ylabel('Speed of Sound, m/s');
 saveas(8,'CO2 Speed of Sound.png')
 %saveas(8,'LaTeX Plots/CO2 Speed of Sound.eps','epsc')

 % Plot the enthalpy vs density
 %??????????????????????????????????????????????????????????????????????????
 figure(9);
 %plot(T,h l(:,1),'b');
 plot(rho l,h l,'??b');hold on; grid on;
 %plot(T,h v(:,1),'r')
 plot(rho v,h v,'r')
 %title('CO2 Enthalpy');
 %legend('Liquid?Helmholtz','Liquid?NIST','Vapor?Helmholtz','Vapor?NIST');
 legend('Liquid?Helmholtz','Vapor?Helmholtz');

 xlabel('Denisty, kg/m3'); ylabel('Enthalpy, kJ/kg');
 saveas(9,'CO2 Enthalpy vs Density.png')
 %saveas(9,'LaTex Plots/CO2 Enthalpy vs Density.eps','epsc')