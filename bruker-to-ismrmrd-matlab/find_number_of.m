function [ number_of_channels, number_of_repetitions , number_of_slices , number_of_echos , number_of_diffusion] = find_number_of( ex )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if(isfield(ex.method,'PVM_EncNReceivers'))
    number_of_channels=ex.method.PVM_EncNReceivers;
else
    number_of_channels=1;
end


if(isfield(ex.method,'PVM_NRepetitions'))
    number_of_repetitions=ex.method.PVM_NRepetitions;
else
    number_of_repetitions=1;
end

if(isfield(ex.method,'PVM_SPackArrNSlices'))
    number_of_slices=ex.method.PVM_SPackArrNSlices;
else
    number_of_slices=1;
end


if(isfield(ex.method,'PVM_NEchoImages'))
    number_of_echos=ex.method.PVM_NEchoImages;
else
    number_of_echos=1;
end


if(isfield(ex.method,'PVM_DwNDiffExp'))
    number_of_diffusion=ex.method.PVM_DwNDiffExp;
else
    number_of_diffusion=1;
end

end

