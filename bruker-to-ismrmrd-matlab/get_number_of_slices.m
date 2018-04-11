function [ number_of_slices ] = get_number_of_slices( ex )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if(isfield(ex.method,'PVM_SPackArrNSlices'))    
    number_of_slices=ex.method.PVM_SPackArrNSlices;
else
    number_of_slices=1;
end

end

