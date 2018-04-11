function [ number_of_averages ] = get_number_of_averages( ex )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if(isfield(ex.method,'PVM_NAverages'))    
    number_of_averages=ex.method.PVM_NAverages;
else
    number_of_averages=1;
end

end

