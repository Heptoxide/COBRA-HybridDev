% Tangential Stress in Cylinder
% using Thick Wall Cylinder Formula

function[stress] = thickCylTan(pressInt, radInt, pressOut, radOut, r)

    stress = (pressInt.*radInt.^2 - pressOut.*radOut.^2 - radInt.^2.*radOut.^2 ...
             *(pressOut-pressInt)./r.^2)./(radOut.^2 - radInt.^2);
             
end