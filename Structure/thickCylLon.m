% Longitudinal Stress in Cylinder
% using Thick Wall Cylinder Formula

function[stress] = thickCylLon(pressInt, radInt, radOut)

    stress = (pressInt.*radInt.^2)./(radOut.^2 - radInt.^2);
             
end