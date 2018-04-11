function [ number_of_repetitions ] = get_number_of_repetitions( ex )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if(isfield(ex.method,'PVM_NRepetitions'))    
    number_of_repetitions=ex.method.PVM_NRepetitions;
else
    number_of_repetitions=1;
end

end

