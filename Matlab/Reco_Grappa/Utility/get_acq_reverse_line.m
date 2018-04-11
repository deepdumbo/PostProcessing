function [ acq_reverse , acq_pas_reverse] = get_acq_reverse_line(meas  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

acq_reverse=  find(  (meas.head.flagIsSet('ACQ_IS_REVERSE')) );

str_msg=sprintf('le nombre de lignes ACQ_IS_REVERSE est  %d \n', size(acq_reverse,2)); disp(str_msg);


acq_pas_reverse=  find(   ~(meas.head.flagIsSet('ACQ_IS_REVERSE')) );

str_msg=sprintf('le nombre de lignes ACQ_IS_PAS_REVERSE est  %d \n', size(acq_pas_reverse,2)); disp(str_msg);


if (size(meas.data,2)~= (size(acq_reverse,2) + size(acq_pas_reverse,2)))
    
    disp('problem number of reverse line');
    return;  
end

end

