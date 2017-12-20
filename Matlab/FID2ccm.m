function [im, W, F] = FID2ccm (mode)

%%%%%        Routine script for reconstruction      %%%%%
%%%%                                                 %%%%
%%%         This script handles the FID reconst.      %%%
%%         to .h5 and a matlab coils combination       %%
%        + phase comparison between Inati & Walsh       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;

%% First add path of the necessary functions
addpath('/home/kylianh/Dev/cute/MatlabSourceCode/CSM/');
addpath('/home/kylianh/Dev/cute/bruker-to-ismrmrd-matlab');
addpath('/home/kylianh/Dev/cute/Imagine/code');

%% FID reconstruction
if(strcmp(mode, 'Reco'))
    disp('--> GENERIC_CREATE_DATASET_FROM_BRUKER');
        generic_create_dataset_from_bruker();
    disp('<-- GENERIC_CREATE_DATASET_FROM_BRUKER');

%% From hdf5 files already reconstructed with generic_create_dataset_from_bruker()
elseif(strcmp(mode, 'h5'))
    disp('--> HDF5_KSPACE_READER');
        data_for_acqp = hdf5_kspace_reader('/Dicom/DIXON/Validation/RecoData/test/',[23, 24, 25]);
    disp('<-- HDF5_KSPACE_READER');
end
    
%% Coils combination - Image reconstruction
disp('--> DIXON_COIL_COMBINE');
    disp('    :: Inati & Walsh reconstruction');
    
        % Number of echoes
        nechoes = 3;
        [im.inati, im.walsh] = Dixon_coil_combine(data_for_acqp, nechoes, mode);
        im.inati = permute(im.inati, [1 2 4 3]);
        im.walsh = permute(im.walsh, [1 2 4 3]);

        figure('name','Inati Reconstruction');
        for s=1:nechoes
            subplot(2,nechoes,s);           imagesc(abs(squeeze(im.inati(:,:,s)))); title(['Echo', num2str(s)]);  axis square;
            subplot(2,nechoes,s+nechoes);   imagesc(angle(squeeze(im.inati(:,:,s)))); axis square;
            colormap(gray);
        end

        figure('name','Walsh Reconstruction');
        for s=1:nechoes
            subplot(2,nechoes,s);           imagesc(abs(squeeze(im.walsh(:,:,s)))); title(['Echo', num2str(s)]);  axis square;
            subplot(2,nechoes,s+nechoes);   imagesc(angle(squeeze(im.walsh(:,:,s)))); axis square;
            colormap(gray);
        end   
disp('<-- DIXON_COIL_COMBINE');

%% 3 points Dixon reconstruction
disp('--> DIXON_3P');
    disp('    :: Dixon processing');
        [W.inati, F.inati, IP.inati, OP.inati] = Dixon_3P(im.inati(:,:,1),im.inati(:,:,2),im.inati(:,:,3), 1);
        [W.walsh, F.walsh, IP.walsh, OP.walsh] = Dixon_3P(im.walsh(:,:,1),im.walsh(:,:,2),im.walsh(:,:,3), 1);

    disp('    :: Water and Fat images');
        figure('name','Dixon :: Inati Reconstruction');
        subplot(121);   imagesc(abs(W.inati)); title('Water Only');   colorbar; axis square;
        subplot(122);   imagesc(abs(F.inati)); title('Fat Only');     colorbar; axis square;
        colormap(gray);

        figure('name','Dixon :: Walsh Reconstruction');
        subplot(121);   imagesc(abs(W.walsh)); title('Water Only');   colorbar; axis square;
        subplot(122);   imagesc(abs(F.walsh)); title('Fat Only');     colorbar; axis square;
        colormap(gray);
disp('<-- DIXON_3P');