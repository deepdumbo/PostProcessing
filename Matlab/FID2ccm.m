function [im, W, F] = FID2ccm (mode)

%%%%%        Routine script for reconstruction      %%%%%
%%%%                                                 %%%%
%%%         This script handles the FID reconst.      %%%
%%         to .h5 and a matlab coils combination       %%
%        + phase comparison between Inati & Walsh       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;

%% First add path of the necessary functions
addpath('CSM/');
addpath('../bruker-to-ismrmrd-matlab');
addpath('Unwrap/quality/');

%% FID reconstruction
if(strcmp(mode, 'Reco'))
    disp('--> GENERIC_CREATE_DATASET_FROM_BRUKER');
        generic_create_dataset_from_bruker();
    disp('<-- GENERIC_CREATE_DATASET_FROM_BRUKER');

%% From hdf5 files already reconstructed with generic_create_dataset_from_bruker()
elseif(strcmp(mode, 'h5'))
    disp('--> HDF5_KSPACE_READER');
        [data_for_acqp, output_tmp] = hdf5_kspace_reader('/Dicom/DIXON/Validation/RecoData/In_Vitro/2D/No_Grappa/20171221/3xFLASH/',[37, 40, 41]);
    disp('<-- HDF5_KSPACE_READER');
end
    
%% Coils combination - Image reconstruction
disp('--> DIXON_COIL_COMBINE');
    disp('    :: Walsh reconstruction');
    
        % Number of echoes
        nechoes = 3;
        im = Dixon_coil_combine(data_for_acqp, nechoes, mode);
        im = permute(im, [1 2 4 3]);
        ismrm_imshow(abs(im),[],[1 nechoes],{'Echo 1' 'Echo 2' 'Echo 3'}, 'Walsh Reconstruction : Magnitude');
        ismrm_imshow(angle(im),[],[1 nechoes],{'Echo 1' 'Echo 2' 'Echo 3'}, 'Walsh Reconstruction : Phase');
disp('<-- DIXON_COIL_COMBINE');

%% 3 points Dixon reconstruction
disp('--> DIXON_3P');
    disp('    :: Dixon processing');
        [W, F] = Dixon_3P(im(:,:,1),im(:,:,2),im(:,:,3), 1);

    disp('    :: Water and Fat images');
          shw(:,:,1) = W;
          shw(:,:,2) = F;
          FDixW = ismrm_imshow(abs(shw),[],[1 2],{'Water' 'Fat'}, 'Dixon :: Walsh Reconstruction');
          clear shw;
disp('<-- DIXON_3P');

%% Save the .jpg
% disp('--> SAVEFIGURE');
%     disp(['    :: ', output_tmp, '_Walsh.jpg']);
%         saveFigure(FRecoW, [output_tmp,'_Walsh.jpg']);
%     disp(['    :: ', output_tmp, '_Walsh_Dixon.jpg']);
%         saveFigure(FDixW, [output_tmp,'_Walsh_Dixon.jpg']);
% disp('<-- SAVEFIGURE');