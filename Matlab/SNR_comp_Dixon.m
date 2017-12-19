function [SNR2, SNR3] = SNR_comp_Dixon (img)

    %%%% This script is to compare SNR in fat only image from Dixon 2 and 3
    %%%% points

    encX = size(img(:,:,1),1);
    encY = size(img(:,:,1),2);

    %% Dixon 2 pts :: Water and Fat images
    W2 = (img(:,:,1) + img(:,:,2)) ./ 2;
    F2 = (img(:,:,1) - img(:,:,2)) ./ 2;

    %% Dixon 3 pts :: Water and Fat images
    [W3, F3, ~, ~] = Dixon_3P(img(:,:,1), img(:,:,2), img(:,:,3), 1);

    %% Define the ROIs
    r = 10;
    id = r+10;
    Noise   = zeros(encX, encY);
    Signal  = zeros(encX, encY);

    Noise((id-r):(id+r),(id-r):(id+r))      = 1;
    Signal((83-r):(83+r),(103-r):(103+r))   = 1;

    %% Get the ROIs data
    nROI2 = createROI(F2, Noise);
    nROI3 = createROI(F3, Noise);
    sROI2 = createROI(F2, Signal);
    sROI3 = createROI(F3, Signal);

    %% SNR calculations
    SNR2 = abs(mean(sROI2(:)) / std(nROI2(:)));
    SNR3 = abs(mean(sROI3(:)) / std(nROI3(:)));

    %% Display images
    r = r - 1;
    Noise((id-r):(id+r),(id-r):(id+r))  = 0;
    Signal((83-r):(83+r),(103-r):(103+r)) = 0;
    dp_ROI = Noise+Signal;

    mnF = mean(abs(F2(:)));
    minF = mnF*0.8;
    maxF = mnF*6;

    figure(1);
    subplot(221);       imagesc(abs(W2+dp_ROI.*max(W2(:)))); title('Water Only (2 pts)');  colorbar;   axis square;
    subplot(222);       imagesc(abs(F2+dp_ROI.*max(F2(:)))); title('Fat Only (2 pts)');    colorbar;   axis square; caxis([minF maxF]);

    mnF = mean(abs(F3(:)));
    minF = mnF*0.8;
    maxF = mnF*6;

    subplot(223);       imagesc(abs(W3+dp_ROI.*max(W3(:)))); title('Water Only (3 pts)');  colorbar;   axis square;
    subplot(224);       imagesc(abs(F3+dp_ROI.*max(F3(:)))); title('Fat Only (3 pts)');    colorbar;   axis square; caxis([minF maxF]);

    colormap(gray);

end

function ROI = createROI(data, mask)

    ROI = data .* mask;
 
    % Crop the 0 from the mask
    ROI(~any(ROI, 2), :)   = []; % rows
    ROI(:, ~any(ROI, 1))   = []; % column

end
