function [ output ] = find_time_difference_between_peaks( input, Tr_ )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
for i=2:size(input);
    
    output(i,1)=round((input(i)-input(i-1))*Tr_);
    output(i,2)=i;
    
end

str_msg=sprintf('mean %f +/- %f \n',  mean(output(:,1))  , std(output(:,1))); disp(str_msg);



end

