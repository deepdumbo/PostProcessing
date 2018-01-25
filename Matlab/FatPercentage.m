% Calculate the fat percentage in a user-defined ROI for Dixon application

function [PercentFat] = FatPercentage(im)
    
    % Works for 2D acquisitions
    assert( size(im,3) > 2, '2-D image must have at least 3 echoes');

    % Study only the abs()
    absim = abs(im);
    
    % In-phase corrected acquisition
    INCOR = ( absim(:,:,1) + absim(:,:,3) ) ./ 2;

    PercentFat.Corrected_Fatmap = ((INCOR - absim(:,:,2))./INCOR);

    [PercentFat.masked_data, PercentFat.mask] = ROImask(PercentFat.Corrected_Fatmap, '% Fat');
    
    PercentFat.min   = min(PercentFat.masked_data(:));
    PercentFat.max   = max(PercentFat.masked_data(:));
    PercentFat.MEAN  = mean(PercentFat.masked_data(:),'omitnan');
    PercentFat.SD    = std(PercentFat.masked_data(:),'omitnan');
    PercentFat.pixel = sum(PercentFat.mask(:) == 1);
    
    disp('-------------------------------------');
    fprintf('  Min/Max: %.2f/%.2f\n  Mean/SD: %.2f/%.2f\n  ROI: %d pixels\n\n  %%Fat: %.2f%%\n',PercentFat.min, PercentFat.max, PercentFat.MEAN, PercentFat.SD, PercentFat.pixel, PercentFat.MEAN*100);
    disp('-------------------------------------');
end


function [masked_data, mask] = ROImask(im, name)
        % Show the image of interest
        image = abs(im);
        RoiFig = figure('Name','Please draw a ROI...' , 'NumberTitle','off'); 
        imagesc(image,[0 1]); title(upper(name)); colormap(gray); colorbar;   zoom off;
        
        % Draw the ROI 
        mask = double(roipoly);
        
        close(RoiFig); 
        
        masked_data = apply_mask(image, mask);
end

function masked_data = apply_mask(image, mask)
        
        masked_data = mask .* image;
        masked_data (mask == 0) = NaN;
end