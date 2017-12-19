function [ output ] = normalize_navigator( input, debut, fin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


min_value=max(input(debut:fin,1));
	
max_value=min(input(debut:fin,1));
	
output=(input-min_value)/(max_value-min_value);

end

