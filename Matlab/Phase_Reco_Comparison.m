%% Comparison of different phase reconstructions :
%       - From Bruker (Shuffle -> no coil combined)
%           -- SWI -> Phase_image
%           -- Walsh reconstruction (matlab code)
%       - From FID (.h5 from generic_matlab)
%           -- Walsh reconstruction (matlab code)
%       - From Gadgetron
%           -- Iterative Inati's method

close all;
clear;

%% On Ex_vivo dataset
    % Initialization
    str_user        = get_PC_name();
    dataset         = 'Validation/RecoData/Ex_Vivo/2D/No_Grappa/20171221/3xFLASH/';
    datafolder      = [19 20 21];
    datatype = 'flash';
%     dataset         = 'Verification/In_Vitro/2D/20180322/';
%     datafolder      = [11:13];
    TE1              = 2.8;
    ES               = 0.35;

    filename        = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/', dataset];
    filename_h5     = ['/home/', str_user, '/Dicom/DIXON/', dataset];
    filename_gad    = ['/home/', str_user, '/Dicom/DIXON/', dataset, num2str(datafolder(1))];
    DEBUG           = 1;
    dimStatus       = '2D';
    
    
if(strcmpi(datatype,'dixon'))
    %% From Bruker (Shuffle -> No coil combined) 
    [Br_Sh, ~] = reco2dseq_Dixon([filename, num2str(datafolder)], 2, 'shuffle');
    Br_Sh = rot90(Br_Sh, 3);
    
    Br_Sh_cplx = complex(Br_Sh(:,:,:,:,:,1), Br_Sh(:,:,:,:,:,2));
    [nX, nY, nZ, nEchoes, nCoils] = size(Br_Sh_cplx);
    
    i = 1:nEchoes;
    TE(i) = TE1 + (i -1)*ES;
    TE = TE';
    
     if(strcmpi(dimStatus,'2D'))
        Br_Sh_cplx = permute(Br_Sh_cplx, [1 2 5 4 3]); % nX nY nCoils nEchoes nZ
     end

    if(~DEBUG)       
        titles = {'Coil 1', 'Coil 2', 'Coil 3', 'Coil 4', 'Coil 5', 'Coil 6', 'Coil 7'};
        ismrm_imshow(abs(Br_Sh_cplx(:,:,1:nCoils,1)),[],[1 nCoils],titles, 'Bruker :: Shuffle - Magnitude');
        ismrm_imshow(angle(Br_Sh_cplx(:,:,1:nCoils,1)),[],[1 nCoils],titles, 'Bruker :: Shuffle - Phase');
    end
    
    % Coil combination : Walsh
    for ne = 1:nEchoes

        % Coil sensibility map estimation
        csm_walsh(:,:,:,ne) = ismrm_estimate_csm_walsh(squeeze(Br_Sh_cplx(:,:,:,ne)));

        % Correct csm to fit the shading profile with a square root sum-of-square channel combination
        csm_walsh(:,:,:,ne) = ismrm_normalize_shading_to_sos(csm_walsh(:,:,:,ne));

        % Computes noise-optimal channel combination maps from  coil sensitivity maps and a noise covariance matrix.
        ccm_walsh(:,:,:,ne) = ismrm_compute_ccm(csm_walsh(:,:,:,ne));

        % Reconstruction of the images
        Br_Walsh(:,:,:,ne) = sum(squeeze(Br_Sh_cplx(:,:,:,ne)) .* ccm_walsh(:,:,:,ne), 3);
    end
    
    Br_Walsh = permute(Br_Walsh,[1 2 4 3]);

    if(DEBUG)
        titles = cellstr(strcat('TE = ',num2str(TE(1:nEchoes,1),'%-.2f'), ' ms'));
        ismrm_imshow(abs(Br_Walsh(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'Bruker :: Walsh - Magnitude');
        ismrm_imshow(angle(Br_Walsh(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'Bruker :: Walsh - Phase');
    end

elseif(strcmpi(datatype,'flash'))
%         %% From Bruker (Shuffle -> No coil combined) 
%         for i = 1:size(datafolder,2)
%             [Br_Sh(:,:,:,:,i), ~]  = reco2dseq([filename, num2str(datafolder(i))], 2, 'shuffle');
%     %         [Br_Swi(:,:,:,:,i), ~] = reco2dseq([filename, num2str(datafolder(i))], 3, 'combined');
%             TE(i) = TE1 + (i-1)*ES;
%         end
% 
%             % Rotate the matrix on the right direction
%             TE = TE';
%             Br_Sh = rot90(Br_Sh, 3);
%     %         Br_Swi = rot90(Br_Swi, 3);
% 
%             % Check 2D or 3D and create complex data
%             if(strcmpi(dimStatus,'2D'))
%                 Br_Sh = permute(Br_Sh, [1 2 3 5 4]);
%                 [nX, nY, nCoils, nEchoes, ~] = size(Br_Sh);
%                 Br_Sh_cplx = complex(Br_Sh(:,:,:,:,1), Br_Sh(:,:,:,:,2));
% 
%     %             Br_Swi = permute(Br_Swi, [1 2 5 4 3]);
%             else
%                 [nX, nY, nZ, nCoils, nEchoes, ~] = size(Br_Sh);
%                 Br_Sh_cplx = complex(Br_Sh(:,:,:,:,1,:), Br_Sh(:,:,:,:,2,:));
%             end
% 
%             % Data are considered 2D (FLASH)
%             if(~DEBUG)       
%                 titles = {'Coil 1', 'Coil 2', 'Coil 3', 'Coil 4', 'Coil 5', 'Coil 6', 'Coil 7'};
%                 ismrm_imshow(abs(Br_Sh_cplx(:,:,1:nCoils,1)),[],[1 nCoils],titles, 'Bruker :: Shuffle - Magnitude');
%                 ismrm_imshow(angle(Br_Sh_cplx(:,:,1:nCoils,1)),[],[1 nCoils],titles, 'Bruker :: Shuffle - Phase');
%             end
% 
%             % Coil combination : Walsh
%             for ne = 1:nEchoes
% 
%                 % Coil sensibility map estimation
%                 csm_walsh(:,:,:,ne) = ismrm_estimate_csm_walsh(squeeze(Br_Sh_cplx(:,:,:,ne)));
% 
%                 % Correct csm to fit the shading profile with a square root sum-of-square channel combination
%                 csm_walsh(:,:,:,ne) = ismrm_normalize_shading_to_sos(csm_walsh(:,:,:,ne));
% 
%                 % Computes noise-optimal channel combination maps from  coil sensitivity maps and a noise covariance matrix.
%                 ccm_walsh(:,:,:,ne) = ismrm_compute_ccm(csm_walsh(:,:,:,ne));
% 
%                 % Reconstruction of the images
%                 Br_Walsh(:,:,:,ne) = sum(squeeze(Br_Sh_cplx(:,:,:,ne)) .* ccm_walsh(:,:,:,ne), 3);
% 
%             end
%                 Br_Walsh = permute(Br_Walsh,[1 2 4 3]);
% 
%             if(DEBUG)
%                 titles = cellstr(strcat('TE = ',num2str(TE(1:nEchoes,1),'%-.2f'), ' ms'));
%     %             ismrm_imshow(Br_Swi(:,:,1:nEchoes),[],[1 nEchoes],titles, 'Bruker :: SWI phase image');
%                 ismrm_imshow(abs(Br_Walsh(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'Bruker :: Walsh - Magnitude');
%                 ismrm_imshow(angle(Br_Walsh(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'Bruker :: Walsh - Phase');
%             end

        %% From FID (.h5 from generic_matlab)
        str_msg = sprintf('\n From FID (Matlab) :'); disp(str_msg);
        for i = 1:size(datafolder,2)
            str_msg = sprintf('\n -- %s', [filename_h5, num2str(datafolder(i))]); disp(str_msg);
            AllData     = demo_flash_grappa_generic_kylian([filename_h5, num2str(datafolder(i))]);
            FID_Walsh(:,:,i)   = rot90(AllData.image.grappa.walsh,1);
        end
            
            % Correct the offset (based on Br_Sh)
            FID_Walsh = circshift(FID_Walsh, 19);
    %           FID_Walsh = circshift(FID_Walsh, 11);
            
            % Data are considered 2D (FLASH)
            if(DEBUG)
                titles = cellstr(strcat('TE = ',num2str(TE(1:nEchoes,1),'%-.2f'), ' ms'));
                ismrm_imshow(abs(FID_Walsh(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'FID :: Walsh - Magnitude');
                ismrm_imshow(angle(FID_Walsh(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'FID :: Walsh - Phase');
            end
        
    %     %% From Gadgetron (Iterative Inati's method)
    %         filename_iter   = [filename_gad, '_', num2str(datafolder(end)), '_img.h5'];
    %         str_msg = sprintf('\n From FID (Gadgetron) :\n -- %s \n', filename_iter); disp(str_msg);
    %         
    %         hinfo           = hdf5info(filename_iter);
    %         res_m           = single(h5read(filename_iter,  hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name));
    %         res_p           = single(h5read(filename_iter,  hinfo.GroupHierarchy.Groups(1).Groups(2).Datasets(2).Name));  
    %         res_p           = (res_p-2048)/2048*pi; 
    %         res_cplx        = res_m.*exp(1i.*res_p);
    %         
    %         Gad_Iter        = rot90(permute(res_cplx,[1 2 5 3 4]), 1);
    % 
    %         clear hinfo res_m res_p res_cplx
    %         
    %         % Correct the offset (based on Br_Sh)
    %         Gad_Iter = circshift(Gad_Iter, 19);
    % %           Gad_Iter = circshift(Gad_Iter, 11);
    %         
    %         % Data are considered 2D (FLASH)
    %         if(DEBUG)
    %             titles = cellstr(strcat('TE = ',num2str(TE(1:nEchoes,1),'%-.2f'), ' ms'));
    %             ismrm_imshow(abs(Gad_Iter(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'Gadgetron :: Inati Iter - Magnitude');
    %             ismrm_imshow(angle(Gad_Iter(:,:,1:nEchoes)),[],[1 nEchoes],titles, 'Gadgetron :: Inati Iter - Phase');
    %         end
end
