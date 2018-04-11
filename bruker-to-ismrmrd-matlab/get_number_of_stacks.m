function [ number_of_stacks ] = get_number_of_stacks( ex )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if(isfield(ex.method,'PVM_NSPacks'))    
    number_of_stacks=ex.method.PVM_NSPacks;
else
    number_of_stacks=1;
end 


end

