function [im, Dixon, TE] = FID2ccm (mode)

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
        echoes  = [30:39];
        [data_for_acqp, output_tmp, TE] = hdf5_kspace_reader('/Dicom/DIXON/Validation/RecoData/In_Vitro/2D/No_Grappa/20180227/',echoes);
        TE = TE';
    disp('<-- HDF5_KSPACE_READER');
end
    
%% Coils combination - Image reconstruction
disp('--> DIXON_COIL_COMBINE');
    disp('    :: Walsh reconstruction');
    
        % Number of echoes
        nshow = 3;
        nechoes = size(TE,1);
        im = Dixon_coil_combine(data_for_acqp, nechoes, mode);
        im = permute(im, [1 2 4 3]);
        titles = cellstr(strcat('TE = ',num2str(TE(1:nshow,1),'%-.2f'), ' ms'));
        FRecoWM = ismrm_imshow(abs(im(:,:,1:nshow)),[],[1 nshow],titles, 'Walsh Reconstruction : Magnitude');
        FRecoWP = ismrm_imshow(angle(im(:,:,1:nshow)),[],[1 nshow],titles, 'Walsh Reconstruction : Phase');
%                 titles = cellstr(strcat('TE = ',num2str(TE(7:nshow,1),'%-.2f'), ' ms'));
%         FRecoWM = ismrm_imshow(abs(im(:,:,7:nshow)),[],[1 3],titles, 'Walsh Reconstruction : Magnitude');
%         FRecoWP = ismrm_imshow(angle(im(:,:,7:nshow)),[],[1 3],titles, 'Walsh Reconstruction : Phase');
disp('<-- DIXON_COIL_COMBINE');

%% 3 points Dixon reconstruction
disp('--> DIXON_3P');
    disp('    :: Dixon processing');
        [W, F] = Dixon_3P(im(:,:,1),im(:,:,2),im(:,:,3),TE);
%         [W, F] = Dixon_3P(im(:,:,7),im(:,:,8),im(:,:,9),TE(7:end));

    disp('    :: Water and Fat images');
          Dixon(:,:,1) = W;
          Dixon(:,:,2) = F;
          figure('Name','Dixon :: Water','NumberTitle','off'); imagesc(abs(W)); axis image; axis off; colormap gray; colorbar; title('Water');
          figure('Name','Dixon :: Fat','NumberTitle','off'); imagesc(abs(F)); axis image; axis off; colormap gray; colorbar; title('Fat');
         % FDixW = ismrm_imshow(abs(Dixon),[],[1 2],{'Water' 'Fat'}, 'Dixon :: Walsh Reconstruction');
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