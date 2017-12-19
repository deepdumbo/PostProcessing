function [ output ] = Filtering( input, mode , Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

output=zeros(size(input));

if (mode==1)

% lowpass Butterworth filter
fNorm = 0.8 / (Fs/2);               %# normalized cutoff frequency
[b,a] = butter(6, fNorm, 'low');  %# 6th order filter
    
elseif (mode==2)
      
fNorm = 40 / (Fs/2);               %# normalized cutoff frequency
[b,a] = butter(6, fNorm, 'low');  %# 6th order filter  
        
elseif (mode==3)
 
fNorm = 10 / (Fs/2);               %# normalized cutoff frequency
[b,a] = butter(6, fNorm, 'low');  %# 6th order filter
    
else
    
end


nCh=size(input,2);

for c=1:1:nCh

 temp  = ButterworthLowpassFilterSixthOrder( input(end:-1:1,c) , a, b);
 temp1  = ButterworthLowpassFilterSixthOrder( temp(end:-1:1) , a, b);
 output(:,c)=temp1;
end


end

