clc; clear all;

 % Creats a plot of C02 temperature vs density for different enthalpies

 % Create a temperature and density vector to calculate CO2 properties at
 %??????????????????????????????????????????????????????????????????????????
 %Tarr = linspace(240,320); % Temperature vector, K (?27.67 to 116.33F)
 T = 250;
 rho = linspace(100,1000); % Density vector, kg/m3
 %rho = 1000;

 % Loop through temperature and density vectors to populate the CO2 properties
 %??????????????????????????????????????????????????????????????????????????
 for i = 1:length(T)
 for j = 1:length(rho)
 Props = CO2Props(T(i),rho(j));
 h(i,j) = Props.h; % Save Enthalpy
 end
 end

 % Create a vector of enthalpies to plot vs density and temperature
 h Known = h(1:15:length(h));
 % Preallocate an array to save the calculated temperature to
 Tsave = zeros(length(h Known),length(rho));

 % Calculate the temperature at the previously defined enthalpies and their
 % corrosponding densities by looping through them and saving temperature

 for iter = 1:length(h Known)
 for k = 1:length(rho)
 Tsave(iter,k) = CO2 rho h 2T(rho(k), h(iter), 300);
 end
 end

 % Plot temperature vs density at each of the enthalpies
 figure(1);
 plot(rho,Tsave(1,:),':b'); hold on; grid on;
 plot(rho,Tsave(2,:),'??r');
 plot(rho,Tsave(3,:),'b');
 plot(rho,Tsave(4,:),'??r');
 plot(rho,Tsave(5,:),'b');
 plot(rho,Tsave(6,:),'??r');
 plot(rho,Tsave(7,:),'b');
 %title('CO2 Temperature vs Density at Specific Enthalpies');
 legend(['h = ' num2str(h Known(1),3) ' kJ/kg'], ...
 ['h = ' num2str(h Known(2),3) ' kJ/kg'], ...
 ['h = ' num2str(h Known(3),3) ' kJ/kg'], ...
 ['h = ' num2str(h Known(4),3) ' kJ/kg'], ...
 ['h = ' num2str(h Known(5),3) ' kJ/kg'], ...
 ['h = ' num2str(h Known(6),3) ' kJ/kg'], ...
 ['h = ' num2str(h Known(7),3) ' kJ/kg'], ...
 'Location','SE');
 xlabel('Density, kg/m3'); ylabel('Temperature, K');
 saveas(1,'CO2 Temp vs Density for Diff Enthalpies.png')