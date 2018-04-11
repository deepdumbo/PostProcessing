function [ number_of_stacks ] = get_number_of_stacks( number_of_slices , slice_acceleration )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if (mod(number_of_slices,slice_acceleration)~=0)
    disp('slice acceleration not compatible with the number of slices');
    return;
else
    number_of_stacks=number_of_slices/slice_acceleration;
end

disp(['number_of_slices : ',num2str( number_of_slices )]); 
disp(['slice_acceleration : ',num2str( slice_acceleration )]); 
disp(['number_of_stacks : ',num2str( number_of_stacks )]); 

end

