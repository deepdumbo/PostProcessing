function [im_walsh] = Dixon_coil_combine( data_for_acqp, nechoes, mode)
% This function aims to compare differents method of coils sensitivity
% correlation
% It therefore reconstructs the image data and merge all the coils images
% We decided to return the Inati solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(mode, 'Reco'))
    % Get the acquisition parameters 
    [ enc_Ny, nCoils, Nx_nechoes ] = size(data_for_acqp);
    enc_Nx = Nx_nechoes / nechoes;

    % Separate the different echoes from the global k-space into nechoes
    % k-space of the same dimensions.
    data_ = zeros(enc_Ny, nCoils, enc_Nx, nechoes);

    for ne = 1:nechoes
        data_(:,:,:,ne) = data_for_acqp(:,:,ne:nechoes:end);
    end

    data_ = permute(data_, [1 3 2 4]);

elseif(strcmp(mode, 'h5'))
    [enc_Ny, enc_Nx, nCoils, ~] = size(data_for_acqp);
    data_ = data_for_acqp;
end

%% Reconstruction - Coils combination 
for ne=1:nechoes
    
    % Data -> image
    data_pre_whitening_remove = squeeze(data_(:,:,:,ne));
    im_pre_whitening_remove = ifft_2D( data_pre_whitening_remove);
    
    % Coils combination
    sum_= zeros(enc_Ny,enc_Nx);
    for c = 1:nCoils
        sum_ = sum_ + abs(im_pre_whitening_remove(:,:,c)).*exp(angle(im_pre_whitening_remove(:,:,c)).*1j);
    end
    
    
    magnitudepre_whitening_remove = sqrt(sum(abs(im_pre_whitening_remove(:,:,:)).^2,3));
    phase_pre_whitening_remove = angle(sum_);
    
    % Coil sensibility map estimation
    csm_walsh(:,:,:,ne) = ismrm_estimate_csm_walsh(ifft_2D(data_pre_whitening_remove));
    %csm_inati(:,:,:,ne) = coil_map_study_2d_Inati( ifft_2D(data_pre_whitening_remove), 5, 3 );
    
    % Correct csm to fit the shading profile with a square root sum-of-square channel combination
    csm_walsh(:,:,:,ne) = ismrm_normalize_shading_to_sos(csm_walsh(:,:,:,ne));
    %csm_inati(:,:,:,ne) = ismrm_normalize_shading_to_sos(csm_inati(:,:,:,ne));
    
    % Computes noise-optimal channel combination maps from  coil sensitivity maps and a noise covariance matrix.
    ccm_walsh(:,:,:,ne) = ismrm_compute_ccm(csm_walsh(:,:,:,ne));
    %ccm_inati(:,:,:,ne) = ismrm_compute_ccm(csm_inati(:,:,:,ne));

    % Reconstruction of the images
    im_walsh(:,:,:,ne) = sum(ifft_2D(data_pre_whitening_remove) .* ccm_walsh(:,:,:,ne), 3);
    %im_inati(:,:,:,ne) = sum(ifft_2D(data_pre_whitening_remove) .* ccm_inati(:,:,:,ne), 3);

%     figure(4)
%     subplot(4,2,1);   imagesc(magnitudepre_whitening_remove); colormap(gray); title('standard'); axis square;
%     subplot(4,2,2);   imagesc(phase_pre_whitening_remove); colormap(gray); axis square;
% 
%     subplot(4,2,3);   imagesc(abs(im_walsh(:,:,:,ne))); colormap(gray);   title('walsh'); axis square;
%     subplot(4,2,4);   imagesc(angle(im_walsh(:,:,:,ne))); colormap(gray); axis square;
% 
%     subplot(3,2,5);   imagesc(abs(im_inati(:,:,:,ne))); colormap(gray);   title('inati'); axis square;
%     subplot(3,2,6);   imagesc(angle(im_inati(:,:,:,ne))); colormap(gray); axis square;
   
end


end

