function [    ] = check_size( vec1, vec2, vec3  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (vec1(1)==vec2(1) &&  vec1(1)==vec3(1))
    
    disp('dimension 1 ok')
else
    
    breaks
end

if (vec1(2)==vec2(2) &&  vec1(2)==vec3(2))
    
    disp('dimension 2 ok')
else
    
    breaks
end

if (vec1(3)==vec2(3) &&  vec1(3)==vec3(3))
    
    disp('dimension 3 ok')
else
    
    breaks
end

end

