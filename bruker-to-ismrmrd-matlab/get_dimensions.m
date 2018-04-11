function [ nX, nY, nZ ] = get_dimensions( ex )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

dimensions=ex.method.PVM_Matrix.*ex.method.PVM_AntiAlias;

nX = dimensions(1);
nY = dimensions(2);

if(size(dimensions,1)>2)
    nZ = dimensions(3);
else
    nZ=1;
end


end

