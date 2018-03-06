function [masked_data, mask] = ROImask(im, name)
        % Show the image of interest
        image = abs(im);
        RoiFig = figure('Name','Please draw a ROI...' , 'NumberTitle','off'); 
        imagesc(image); title(upper(name)); colormap(gray); colorbar;   zoom off;
        
        % Draw the ROI 
        mask = double(roipoly);
        
        close(RoiFig); 
        
        masked_data = apply_mask(image, mask);
end

function masked_data = apply_mask(image, mask)
        
        masked_data = mask .* image;
        masked_data (mask == 0) = NaN;
end