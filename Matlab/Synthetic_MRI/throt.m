function Rth=throt(phi,theta)

Rz = zrot(-theta);
Rx = xrot(phi);
Rth = Rz\Rx*Rz; % Inv(Rz)*Rx * Rz


