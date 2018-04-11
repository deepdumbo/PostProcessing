function [ number_of_lines_in_acs , number_of_lines_in_sb ] = get_number_of_lines_in_acs_and_sb( meas )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

acq_parallel = find((meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));

rep=1;

acq_single_band = find( (meas.head.idx.repetition==(rep-1)) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))    );

number_of_lines_in_acs=size(acq_parallel,2);
number_of_lines_in_sb=size(acq_single_band,2);

disp(['number_of_lines_in_acs : ',num2str( number_of_lines_in_acs )]); 
disp(['number_of_lines_in_sb : ',num2str( number_of_lines_in_sb )]); 

end

