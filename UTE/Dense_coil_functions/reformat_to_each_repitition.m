function [ output ] = reformat_to_each_repitition( input  , repetition)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tempo=reshape(input, [size(input,1)/repetition , repetition, size(input,2) ]);
output=permute(tempo, [1 3 2]);
end

