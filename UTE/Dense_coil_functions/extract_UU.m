function [ UU ] = extract_UU( matrixG , debug)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


nbPointsZero=size(matrixG,3);

UU=zeros(size(matrixG,1),size(matrixG,2),nbPointsZero);

for npZ=1:nbPointsZero
    
    A=matrixG(:,:,npZ);
    
    size(A);
%     [U,D] = eig(A);
    
    [U,S,V] = svd(A,'econ');
    
    % http://stackoverflow.com/questions/13704384/does-matlab-eig-always-returns-sorted-values
    % [U,D] = eigs(A,size(A,1));
           
     if (strcmp(debug,'Y'))
    
    figure(6)
    subplot(2,3,(npZ-1)*3+1);imagesc(matrixG(:,:,npZ), [0 1]); colormap(gray);
    subplot(2,3,(npZ-1)*3+2);imagesc(U, [0 1]); colormap(gray);
    subplot(2,3,(npZ-1)*3+3);imagesc(D, [0 1]); colormap(gray);
    
     end
    
    UU(:,:,npZ)=U;
    
end


end

