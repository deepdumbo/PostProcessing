function [ nX, nY, nZ ] = get_dimensions_acq( ex )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


dimensions=ex.method.PVM_EncMatrix;

nX = dimensions(1);
nY = dimensions(2);

if(size(dimensions,1)>2)
    nZ = dimensions(3);
else
    nZ=1;
end


end

