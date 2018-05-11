%%% This function calculates the SNR with ROIs
%%% It takes the same ROI in signal as ROIs for T2* mapping

function [SNR] = SNR_ROI(name, im, mask_noise, mask_signal)
    
    %% Création du masque - Mise en forme des data
    switch nargin
        case 2
            [signal.data, signal.mask]   = ROImask(im, name);
            [noise.data , noise.mask ]   = ROImask(im, 'Noise Area');
        case 3
            noise.data                  = apply_mask(abs(im), mask_noise);
            noise.mask                  = mask_noise;
            [signal.data , signal.mask ]  = ROImask(im, name);
        case 4
            signal.data                 = apply_mask(abs(im), mask_signal);
            signal.mask                 = mask_signal;
            noise.data                  = apply_mask(abs(im), mask_noise);
            noise.mask                  = mask_noise;
        otherwise
            error('Valid number of arguments: 2-4');
    end
    
    %% SNR calculation
        % Mean of signal
        sigmean  = mean(signal.data(:), 'omitnan');
        
        % Std of noise
        noisestd = std(noise.data(:), 'omitnan');
        
        % SNR
        SNR.val         = sigmean ./ noisestd;
        SNR.signal      = signal;
        SNR.noise       = noise ; 

end

function masked_data = apply_mask(image, mask)
        
        masked_data = mask .* image;
        masked_data (mask == 0) = NaN;
end