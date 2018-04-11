function [ index ] = map_interleaved_acquisition( number_of_slices, no_reordering )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(no_reordering,'Y'))
        
    index=0:number_of_slices-1;
    
else
    
    if mod(number_of_slices,2)==1
        index(1)=0;
    else
        index(1)=1;
    end
    
    for i=2:number_of_slices
        index(i)= index(i-1)+2;
        if index(i)>=number_of_slices
            if mod(number_of_slices,2)==0
                index(i)=0;
            else
                index(i)=1;
            end
        end
    end
    
end

% disp('ordre d acquisition des coupes');
% index

end

