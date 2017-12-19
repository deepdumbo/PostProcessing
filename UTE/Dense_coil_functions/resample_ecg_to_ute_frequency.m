function [ output ] = resample_ecg_to_ute_frequency( input , size_limit_output)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

ratio=size_limit_output/size(input,1);

for j=1:1:size(input,2)
    
v=input(:,j);
newNum=size(v,1)*ratio
X =  linspace(0,1,numel(v));
Xi = linspace(0,1,newNum);
iv = interp1(X, v, Xi, 'linear');
output(:,j)=iv';

end

output=single(output);

end

