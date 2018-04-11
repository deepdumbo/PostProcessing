function [ MB_factor , Blipped_CAIPI] = read_wip_user_parameter( temp_struct_user_double )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for ll=1:size(temp_struct_user_double,2)
    
    if (strcmp(temp_struct_user_double(ll).name,'MB_factor'))
        MB_factor=temp_struct_user_double(ll).value;
        str_msg=sprintf('%s %f', temp_struct_user_double(ll).name , temp_struct_user_double(ll).value); disp(str_msg);
    end
    if (strcmp(temp_struct_user_double(ll).name,'Blipped_CAIPI'))
        Blipped_CAIPI=temp_struct_user_double(ll).value;
        str_msg=sprintf('%s %f', temp_struct_user_double(ll).name , temp_struct_user_double(ll).value); disp(str_msg);
    end
    
   
end




end

