function [ readout, E1, E2 ] = get_encoding_size( ex , nZ )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here


readout=ex.method.PVM_EncMatrix(1);

E1=ex.method.PVM_EncMatrix(2);

if(nZ>1)
    E2=ex.method.PVM_EncMatrix(3);
else
    E2=1;
end

end

