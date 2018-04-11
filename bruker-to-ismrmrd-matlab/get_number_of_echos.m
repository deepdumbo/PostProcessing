function [ number_of_echos ] = get_number_of_echos( ex )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if(isfield(ex.method,'PVM_NEchoImages'))
    number_of_echos=ex.method.PVM_NEchoImages;
else
    number_of_echos=1;
end

end

