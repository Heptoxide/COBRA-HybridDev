function result = aerodrag(u, v, w, h)

%	This function gives a shape to the drag curve.
%  This is a first order function at best.

%	VARIABLES and UNITS
%
%	conventions:	Large arrays are in all caps
%						Vectors and other variables are in mixed case, as appropriate
%
%	h					:Height above sea level, in m.
%	u, v, w		:Component velocities of the vehicle (wrt Body coordinates), in m/s.
%	V					:Velocity vector for the vehicle, in m/s.
%	cs_area			:Frontal cross-sectional area of the vehicle, in m2.
%	Cd					:Coefficient of drag, unitless.
%	M					:Mach number the vehicle is at, unitless.
%	gamma				:Coefficient of specific heats for air, unitless.
%	R					:Gas constant for air, in J/kg/K.



global data;
global index;


%	CONSTANTS
if index == 1
   diameter = 	data(1, 47);
	gamma = 		data(1, 97);					% unitless
   R = 			data(1, 98);					% J/kg/K
else
	diameter = 	data(1,       47);
	gamma = 		data(index-1, 97);					% unitless
   R = 			data(index-1, 98);					% J/kg/K
end


cs_area = pi * (diameter/2)^2;				% m2



%	MAIN PROGRAM

V = [u v w];

a = sqrt(gamma * R * temperature(h));
M = norm(V) / a;

Cd = 0.15;
%if (norm(V) <= 100)
%   Cd = .15;
   
%elseif (norm(V) <= 200)
%   Cd = .35;
   
%elseif (norm(V) <= 250)
%   Cd = 1.10;
   
%elseif (norm(V) <= 300)
%   Cd = 1.25;
   
%elseif (norm(V) <= 500)
%   Cd = 1.15;
   
%elseif (norm(V) <= 800)
%   Cd = .95;
   
%else 
%   Cd = .80;
   
%end


Fdx = -1 * sign(u) * Cd * 1/2 * density(h) * u^2 * cs_area;
Fdy = -1 * sign(v) * Cd * 1/2 * density(h) * v^2 * cs_area;
Fdz = -1 * sign(w) * Cd * 1/2 * density(h) * w^2 * cs_area;

Mdx = 0;
Mdy = 0;
Mdz = 0;



result = [Fdx Fdy Fdz Mdx Mdy Mdz];


% Record data
data(index, 97) = gamma;
data(index, 98) = R;
data(index,105) = Cd;
data(index,106) = Fdx;
data(index,107) = Fdy;
data(index,108) = Fdz;
data(index,109) = Mdx;
data(index,110) = Mdy;
data(index,111) = Mdz;

