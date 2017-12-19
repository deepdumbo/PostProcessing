function [  ] = check_if_zero_problem( ex )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

[ nX, nY, nZ ] = get_dimensions( ex );

[ readout, E1, E2 ] = get_encoding_size( ex , nZ )

number_of_repetitions=ex.method.PVM_NRepetitions;
number_of_slices=ex.method.PVM_SPackArrNSlices;
number_of_averages=ex.method.PVM_NAverages;
number_of_echos=ex.method.PVM_NEchoImages;
number_of_stacks=ex.method.PVM_NSPacks;

number_of_channels=ex.method.PVM_EncNReceivers;

pre_calculated_size=readout*number_of_channels*E1*E2*number_of_echos*number_of_repetitions*number_of_slices*number_of_averages*number_of_stacks;



if (mod(size(ex.fid,1)/readout,number_of_channels)==0)
    
    disp('pas de zero en trop');
    remove_the_supplementary_zero=0;
    fprintf('en effet size fid %d  >  %d  (pre-calculated size) \n', size(ex.fid,1), pre_calculated_size );
    
else
    
    disp('il y a des zero en trop');     
    
    if (size(ex.fid,1)> pre_calculated_size)
        
        fprintf('en effet size fid %d  >  %d  (pre-calculated size) \n', size(ex.fid,1), pre_calculated_size );
        remove_the_supplementary_zero=1;
    else
        
        disp('il y a en plus un autre probleme à résoudre');
        fprintf('en effet size fid %d  <=  %d  (pre-calculated size) \n', size(ex.fid,1), pre_calculated_size );
        remove_the_supplementary_zero=1;
    end
    
end
        



end

