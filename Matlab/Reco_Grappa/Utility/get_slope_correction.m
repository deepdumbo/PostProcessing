function [ corrpos_, corrneg_ , tvec ] = get_slope_correction( input_navigator )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% ifft must be done this function
% dimension
% 1 readout
% 2 navigator EPI 1 & 2 & 3
% 3 channels

readout=size(input_navigator,1);    
corrpos_=zeros( readout,1);
corrneg_=zeros( readout,1);
 
[ slope, intercept, x ] = fit_slope( input_navigator );
      
%% application des coefficients
tvec = slope*x + intercept;
str_msg=sprintf('slope %f %f %f \n',  slope, intercept, tvec(1)); disp(str_msg);
    
corrpos_ = exp(complex(zeros( readout, 1), -1*tvec));
corrneg_ = exp(complex(zeros( readout, 1), +1*tvec));

end

