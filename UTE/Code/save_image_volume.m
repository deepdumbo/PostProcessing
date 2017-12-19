function [  ] = save_image_volume(volume, folder_save,  number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  filter )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


str_cardiac=sprintf('%d',number_of_output_cardiac_phases);
str_resp=sprintf('%d',number_of_output_respiratory_phases);

filename_save = [folder_save ,  'image_resp_', str_resp  ,'_cardiac_' ,str_cardiac , '_filter_', filter , '.mat' ];

disp(filename_save);

save(filename_save,'volume', '-v7.3' );


end

