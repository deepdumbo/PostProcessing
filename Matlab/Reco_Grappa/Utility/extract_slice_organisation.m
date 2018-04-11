function [order_of_acquisition_mb, order_of_acquisition_sb, MapSliceSMS, vec_MapSliceSMS]=extract_slice_organisation(number_of_slices,number_of_stacks, MB_factor,  no_reordering)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reordering
% Assuming interleaved SMS acquisition, define the acquired slices for each
% multislice band
MapSliceSMS={};
vec_MapSliceSMS=[];

if (number_of_stacks==1)
    
    str_msg=sprintf('number of stack %d ', number_of_stacks); disp(str_msg);
    MapSliceSMS{1}= map_interleaved_acquisition(MB_factor, no_reordering);
    order_of_acquisition_mb=map_nterleaveACQ_MB(number_of_stacks, no_reordering);
    order_of_acquisition_sb=map_interleaved_acquisition(number_of_slices, no_reordering);    
    vec_MapSliceSMS=MapSliceSMS{1};
    
else
    
    str_msg=sprintf('number of stack %d ', number_of_stacks); disp(str_msg);
    
    order_of_acquisition_mb=map_interleaved_acquisition(number_of_stacks, no_reordering);
    order_of_acquisition_sb=map_interleaved_acquisition(number_of_slices, no_reordering);
    
    for i=1:number_of_stacks
        
        str_msg=sprintf('stack %d  qui correspond à la coupe geometrique %d ', i, order_of_acquisition_mb(i)); disp(str_msg);
        
        MapSliceSMS{i}= order_of_acquisition_mb(i):number_of_stacks:number_of_slices-1;
        
        for j=1:MB_factor
            
            str_msg=sprintf('%d  ', MapSliceSMS{i}(j)); disp(str_msg);
            vec_MapSliceSMS=[vec_MapSliceSMS , MapSliceSMS{i}(j)+1];
        end
        
        %         str_msg=sprintf('\n'); disp(str_msg);
    end
    
   % order_of_acquisition_mb=mapInterleaveACQ_MB(number_of_stacks, no_reordering);
    
end
