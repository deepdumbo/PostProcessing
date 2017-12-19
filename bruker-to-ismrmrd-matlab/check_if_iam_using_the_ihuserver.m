function [ output ] = check_if_iam_using_the_ihuserver( input )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

output=0;

if (strcmp(input,'valeryo'))
    output=1;
elseif (strcmp(input,'juliem'))
    output=1;
elseif (strcmp(input,'pierreb'))
    output=1;
else
    output=0;
    
    if(output==1)
       disp('je suis sur le serveur'); 
    end
    
end





