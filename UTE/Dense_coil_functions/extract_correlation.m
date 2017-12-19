function [ matrixC, matrixG,matrixS ] = extract_correlation( input , thresholdC,  debug )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nCh=size(input,2);

matrixC=zeros(nCh,nCh);
matrixG=zeros(nCh,nCh);
matrixS=zeros(nCh,nCh);



for c1=1:nCh
    for c2=1:1:nCh
        
        [ r ] = correlation2( squeeze(input(:,c1)) , squeeze(input(:,c2) ));
        
        %             r= xcov(squeeze(ReducerawDataCh(npZ,c1,:)) , squeeze(ReducerawDataCh(npZ,c2,:)),0,'coef');
        
        matrixC(c1,c2)=r;
        %matrixC(c2,c1,npZ)=matrixC(c1,c2,npZ);
        
        if (abs(matrixC(c1,c2))>thresholdC)
            matrixG(c1,c2)=1;
        else
            matrixG(c1,c2)=0;
        end;
        
        if (matrixC(c1,c2)>0)
            matrixS(c1,c2)=1;
        else
            matrixS(c1,c2)=-1;
            
        end
        
        
    end;
end;

if (strcmp(debug,'Y'))
    
    figure(2)
    
    subplot(2,3,(npZ-1)*3+1);imagesc(matrixC(:,:), [0 1]); colormap(gray);
    subplot(2,3,(npZ-1)*3+2);imagesc(matrixG(:,:), [0 1]); colormap(gray);
    subplot(2,3,(npZ-1)*3+3);imagesc(matrixS(:,:), [0 1]); colormap(gray);
    
end




end

