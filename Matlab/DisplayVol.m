function DisplayVol(W, F, IP, OP)
% This function displays volumes of DIXON_3P processing
% in an animation

    aW = abs(W);
    aF = abs(F);
    aIP = abs(IP);
    aOP = abs(OP);

    % All the images have the same size (come from the same dataset)
    [~, ~, dimz] = size(aW);

    close(figure(1))
    figure(1)

    for z = 1:dimz
        
        subplot(221);
        imagesc(aW(:,:,z)); colormap(gray);
        
        subplot(222);
        imagesc(aF(:,:,z));
        
        subplot(223);
        imagesc(aIP(:,:,z));
        
        subplot(224);
        imagesc(aOP(:,:,z));
        
        pause(0.05);
    end


end

