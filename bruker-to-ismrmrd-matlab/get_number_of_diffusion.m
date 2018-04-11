function [ number_of_diffusion ] = get_number_of_diffusion( ex )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


if(isfield(ex.method,'PVM_DwNDiffExp'))
    number_of_diffusion=ex.method.PVM_DwNDiffExp;
else
    number_of_diffusion=1;
end

end

