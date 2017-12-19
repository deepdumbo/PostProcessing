function [ output ] = cut_ecg_at_the_end_of_the_sequence( input )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

compteur=0;

for i=1:1:size(input,1)-6
    
    if(input(i,2)==0 && input(i+1,2)==0 && input(i+2,2)==0 && input(i+3,2)==0 && input(i+4,2)==0)   
    
    compteur=compteur+1;
    if(compteur==1)       
      cut_ecg=i;
        end
    end
        
end

output=input(1:cut_ecg,:);



end

