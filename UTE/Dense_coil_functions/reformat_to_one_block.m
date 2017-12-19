function [ output ] = reformat_to_one_block( input )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



tempo=permute(input,[2 1 3]);

tempo1=reshape(tempo, [size(tempo,1), size(tempo,2) *size(tempo,3) ]);

output=permute(tempo1,[2 1 ]);

end

