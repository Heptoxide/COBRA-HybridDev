% ++_____________________________________________________________________________++
%  This MATLAB function will perform a Runge-Kutta fourth order integration given
%  a specified function name and the function arguments.  The function called must
%  return the partial derivatives (PDR) of the desired equations to integrate.
%  The current version of this routine can handle up to 10 function arguments.
%  Note that the PDRs must be returned as a column array.
%
%  Inputs:
%	PDRFuncName	string		Name of the PDR function
%	TimeArr		array(1,*)	Row array of time points
%	InitCond	array(*,1)      Column array of initial conditions
%					  (same number of rows as PDRs)
%	PDRFuncArg*	any		Required PDR function arguments
%
%  Outputs:
%	TimeArr		array(1,*)	Same as input
%	Y		array(m,n)	Array containing the integrated
%					  values over the time interval
%					     m -> number of PDR equations
%					     n -> number of time steps
%
%					Auth/date: P. Shriver 6/16/99
% ++_____________________________________________________________________________++

function [TimeArr,Y]=Rk4(PDRFuncName,TimeArr,InitCond,PDRFuncArg1,PDRFuncArg2, ...
			 PDRFuncArg3,PDRFuncArg4,PDRFuncArg5,PDRFuncArg6, ...
			 PDRFuncArg7,PDRFuncArg8,PDRFuncArg9,PDRFuncArg10)


%  create the string for the function arguments (to be used in the function call)
	PDRFuncArgs = [];
	for I = 1:nargin-3,
	   PDRFuncArgs=[PDRFuncArgs,' ,PDRFuncArg',int2str(I)];
	end


%  determine the time steps, initialize the Yi array to the initial conditions, and
%  determine the time interval widths
	SizeTimeArr=size(TimeArr);
	NumSteps=SizeTimeArr(2);
	Yi=InitCond';
	Time=TimeArr(1);
	DeltaTime=TimeArr(2)-TimeArr(1);


%  numerically integrate over the time steps, via the PDRs from the function
	for I = 1:NumSteps,
		Y0 = Yi;      % set Y0 to the current time step values
		K0 = eval([PDRFuncName '(Time,Yi' PDRFuncArgs ')']);

		Time = Time + 0.5 * DeltaTime;  % step forward 1/2 time step
		Yi = Y0 + 0.5 * DeltaTime * K0;
		K1 = eval([PDRFuncName '(Time,Yi' PDRFuncArgs ')']);

		% K2 is determined @ same time as K1
		Yi = Y0 + 0.5 * DeltaTime * K1;
		K2 = eval([PDRFuncName '(Time,Yi' PDRFuncArgs ')']);  

		Time = Time + 0.5 * DeltaTime;  % K3 determined @ full time step
		Yi = Y0 + DeltaTime * K2;
		K3 = eval([PDRFuncName '(Time,Yi' PDRFuncArgs ')']);

		% calculate the values at the next full time step
		Yi = Y0 + 1./6. * (K0 + 2.*(K1+K2)+K3) * DeltaTime;
	
		Y(I,:) = Yi';
	end
   Time = Time';
