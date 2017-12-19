function [ img_comb ] = coil_combine( rawDataECG , FT, N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    img_comb=single(zeros(N,N,N));   

    nCh=size(rawDataECG,3);

%% Coil combine
    matricetemp = reshape(rawDataECG,[],nCh);
    for c=1:nCh
        img_sens(:,:,:,c) = FT'*matricetemp(:,c);
        img_comb = img_comb + pow2(abs(img_sens(:,:,:,c)),1);        
    end
    
    img_comb=abs(img_comb);    
   
    
end

