function [im_3D] = Dixon_coil_combine_3D( data_for_acqp, nechoes, mode)
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
        data_ = zeros(enc_Nx, nCoils, enc_Ny, nechoes);

        for ne = 1:nechoes
            data_(:,:,:,ne) = data_for_acqp(:,:,ne:nechoes:end);
        end

        data_ = permute(data_, [1 3 2 4]);

    elseif(strcmp(mode, 'h5'))
        [enc_Nx, enc_Ny, enc_Nz, nCoils, ~] = size(data_for_acqp);
        data_ = data_for_acqp;
        data_ = permute(data_,[1 2 4 5 3]);
    end
    
    hne     = waitbar(0, 'Echo : ', 'Name','Coil combination');
    hz      = waitbar(0, 'Slice : ', 'Name','Coil combination');
    
    pos_w1  =get(hne,'position');
    pos_w2  =[pos_w1(1) pos_w1(2)+pos_w1(4) pos_w1(3) pos_w1(4)];
    set(hz,'position',pos_w2,'doublebuffer','on')
    
    for z=1:enc_Nz
        
        % Update waitbar
            step = z / enc_Nz;
            waitbar(step, hz , ['Slice : ',num2str(z),'/',num2str(enc_Nz)]);
            
        tmp = squeeze(data_(:,:,:,:,z));
        %% Reconstruction - Coils combination 
        for ne=1:nechoes
            
            % Update waitbar
            stepne = ne / nechoes;
            waitbar(stepne, hne , ['Echo : ',num2str(ne),'/',num2str(nechoes)]);
            
            % Data -> image
            data_pre_whitening_remove = squeeze(tmp(:,:,:,ne));
            im_pre_whitening_remove = ifft_2D( data_pre_whitening_remove);

            % Coils combination
            sum_= zeros(enc_Nx,enc_Ny);
            for c = 1:nCoils
                sum_ = sum_ + abs(im_pre_whitening_remove(:,:,c)).*exp(angle(im_pre_whitening_remove(:,:,c)).*1j);
            end

            % Coil sensibility map estimation
            csm_walsh(:,:,:,ne) = ismrm_estimate_csm_walsh(ifft_2D(data_pre_whitening_remove));

            % Correct csm to fit the shading profile with a square root sum-of-square channel combination
            csm_walsh(:,:,:,ne) = ismrm_normalize_shading_to_sos(csm_walsh(:,:,:,ne));

            % Computes noise-optimal channel combination maps from  coil sensitivity maps and a noise covariance matrix.
            ccm_walsh(:,:,:,ne) = ismrm_compute_ccm(csm_walsh(:,:,:,ne));

            % Reconstruction of the images
            im_walsh(:,:,:,ne) = sum(ifft_2D(data_pre_whitening_remove) .* ccm_walsh(:,:,:,ne), 3);

        end
            im_3D(:,:,:,z) = permute(im_walsh, [1 2 4 3]); 
    end
    
    close(hz);
    close(hne);
end
