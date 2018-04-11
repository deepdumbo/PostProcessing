function [ number_of_lines_in_image ,  number_of_lines_in_navigator , number_of_lines_in_image_reverse ] = get_number_number_of_lines_nav_and_image(meas , is_wip_sms)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


contrast=1;

for rep=1:1:1
    
    for s=1:1:1
        
        acq_test = find(  (meas.head.idx.contrast==(contrast-1)) ...
            & (meas.head.idx.slice==(s-1)) ...
            & (meas.head.idx.repetition==(rep-1))...
            & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));
        
        acq_test   = find(  (meas.head.idx.contrast==(contrast-1)) ...
            & (meas.head.idx.repetition==(rep-1)) ...
            & (meas.head.idx.slice==(s-1))...
            & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
            & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))    )  ;
        
        number_of_lines_in_navigator=size(acq_test,2);
        
        
        acq_test   = find(  (meas.head.idx.contrast==(contrast-1)) ...
            & (meas.head.idx.repetition==(rep-1)) ...
            & (meas.head.idx.slice==(s-1))...
            & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
            & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))    )  ;
        
        number_of_lines_in_image=size(acq_test,2);
        
        
        acq_test   = find(  (meas.head.idx.contrast==(contrast-1)) ...
            & (meas.head.idx.repetition==(rep-1)) ...
            & (meas.head.idx.slice==(s-1))...
            & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
            & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  ...
            & (meas.head.flagIsSet('ACQ_IS_REVERSE'))   );
        
        number_of_lines_in_image_reverse=size(acq_test,2);
        
    end
    
end

disp(['number_of_lines_in_navigator : ',num2str( number_of_lines_in_navigator )]); 
disp(['number_of_lines_in_image : ',num2str( number_of_lines_in_image )]); 
disp(['number_of_lines_in_image_reverse : ',num2str( number_of_lines_in_image_reverse )]); 

if (strcmp(is_wip_sms,'Y'))
    
    number_of_lines_in_navigator=number_of_lines_in_navigator/2;
    number_of_lines_in_image=number_of_lines_in_image/2;
    number_of_lines_in_image_reverse= number_of_lines_in_image_reverse/2;
end


disp(['number_of_lines_in_navigator : ',num2str( number_of_lines_in_navigator )]); 
disp(['number_of_lines_in_image : ',num2str( number_of_lines_in_image )]); 
disp(['number_of_lines_in_image_reverse : ',num2str( number_of_lines_in_image_reverse )]); 

end

