function [im, Dixon, TE] = FID2ccm_3D (mode)

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
        echoes  = [14 20 14];
        
        if(size(echoes,2) > 1)
            [data_for_acqp, output_tmp, TE] = hdf5_kspace_reader_3D('/Dicom/DIXON/Sheep_FAT/',echoes);
        else
            [data_for_acqp, output_tmp, TE] = hdf5_kspace_reader_multi_echoes('/Dicom/DIXON/Validation/RecoData/Ex_Vivo/2D/No_Grappa/20180216/',echoes);
        end
        
        TE = TE';
    disp('<-- HDF5_KSPACE_READER');
end
    
%% Coils combination - Image reconstruction
disp('--> DIXON_COIL_COMBINE');
    disp('    :: Walsh reconstruction');
    
        % Number of echoes
        nshow = 3;
        nechoes = size(TE,1);
        
        im = Dixon_coil_combine_3D(data_for_acqp, nechoes, mode);
        
        size(im)
        titles = cellstr(strcat('TE = ',num2str(TE(1:nshow,1),'%-.2f'), ' ms'));
        FRecoWM = ismrm_imshow(abs(im(:,:,1:nshow)),[],[1 nshow],titles, 'Walsh Reconstruction : Magnitude');
        FRecoWP = ismrm_imshow(angle(im(:,:,1:nshow)),[],[1 nshow],titles, 'Walsh Reconstruction : Phase');
        
        im = permute(im, [1 2 4 3]);
        Dixon = 0;
disp('<-- DIXON_COIL_COMBINE');

%% 3 points Dixon reconstruction
disp('--> DIXON_3P');
    disp('    :: Dixon processing');
        [W, F] = Dixon_3P(im(:,:,1),im(:,:,2),im(:,:,3),TE);

    disp('    :: Water and Fat images');
          Dixon(:,:,1) = W;
          Dixon(:,:,2) = F;
          FDixW = ismrm_imshow(abs(Dixon),[],[1 2],{'Water' 'Fat'}, 'Dixon :: Walsh Reconstruction');
disp('<-- DIXON_3P');

% %% Save the .jpg
% disp('--> SAVEFIGURE');
%     disp(['    :: ', output_tmp, '_Mag_Walsh.jpg']);
%         saveFigure(FRecoWM, [output_tmp,'_Mag_Walsh.jpg']);
%     disp(['    :: ', output_tmp, '_Phase_Walsh.jpg']);
%         saveFigure(FRecoWP, [output_tmp,'_Phase_Walsh.jpg']);
%     disp(['    :: ', output_tmp, '_Walsh_Dixon.jpg']);
%         saveFigure(FDixW, [output_tmp,'_Walsh_Dixon.jpg']);
% disp('<-- SAVEFIGURE');