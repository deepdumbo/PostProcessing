%%%%                        Script to get a csm simulation
%%%%
%%%% Protocole : 
%%%%        1. Generate a .h5 from /usr/local/bin/ismrmrd_generate_cartesian_shepp_logan
%%%%           with matrix dim, ncoils, noise desired (matrix dim should
%%%%           match encX, encY from Generate_Synthetic_FLASH_image)
%%%%        2. Extract the csm from the .h5 file
%%%%        3. Simulate FLASH image with ncoils with/without noise
%%%%        4. Go to ../../Matlab/FID2ccm and apply the DIXON

close all;
clear all;

disp('--> GENERATE_NOISY_MULTICOILS_FLASH_SIMULATION');
    disp('    :: Dixon Fat/Water application');
        addpath('../');
        addpath('../CSM');
        addpath('../Unwrap/quality');
        addpath('../../ismrm_sunrise_matlab/');

        %% Load session data
        [~,id]= system('whoami');
        str_user= id(1:end-1);       
        clear id
        
        %% Parameters 
        acc_factor = 1;
        noise = 0.05;
        
        %% Load the filename and data
        filename = 'example.h5';

        if exist(filename, 'file')
            hinfo = hdf5info(filename);   
            smaps = h5read(filename,  hinfo.GroupHierarchy.Groups(1).Datasets(2).Name);  
            smaps = complex(smaps.real, smaps.imag);
%             smaps=circshift(smaps,-2,3);
            ncoils = size(smaps,3);
        else
            error(['File ' filename ' does not exist.  Please generate it.'])
        end  
        
        % Generate no noise synthetic phantom
        Generate_Synthetic_FLASH_images;

        % Check if simulation is 2D
        if(size(synFLASH,3) == 1)
            synFLASH = permute(synFLASH, [1 2 4 3]);
        end        

        %% Use of ismrm_hansen_generate_data.m by Michael S. Hansen but softly modified
    disp('    :: Noise application & 7 channels simulation');

%         load('noise_covariances.mat');
%         Rn_normal_7 = Rn_normal_8(1:ncoils,1:ncoils);
%         clear Rn_normal_8
        
        Rn_normal_7 = eye(7,7); 
        L_normal_7 = chol(Rn_normal_7,'lower');

        for i=1:ne
            img = synFLASH(:,:,i);
            noise_level = noise*max(img(:));

            %Sample Cartesian Data
            noise_white = noise_level*complex(randn(size(img,1),size(img,2),size(smaps,3)),randn(size(img,1),size(img,2),size(smaps,3)));
            noise_color = reshape( permute( L_normal_7 * permute( reshape( noise_white, numel(noise_white)/ncoils,ncoils ),[2 1] ),[2 1] ), size(noise_white) );

            [data, sp] = ismrm_sample_data(img, smaps, acc_factor, 0);
            synData(:,:,:,i) = data + noise_color .* repmat(sp > 0,[1 1 ncoils]);
        end
        
        titles = cellstr(strcat('Coil',num2str([1:size(synData,3)]','%d')));
        ismrm_imshow(abs(ifft_2D(squeeze(synData(:,:,:,1)))),[],[1 size(synData,3)],titles,'Synthetic MRI :: synData coils');

        %% Coils combination - Image reconstruction
        disp('--> DIXON_COIL_COMBINE');
            disp('    :: Walsh reconstruction');

                % Number of echoes
                nechoes = size(TE,1);
                im = Dixon_coil_combine(synData, nechoes, 'h5');
                im = permute(im, [1 2 4 3]);
                
                titles = cellstr(strcat('TE = ',num2str(TE,'%-.2f'), ' ms'));
                ismrm_imshow(abs(im),[],[1 size(im,3)],titles,'Synthetic MRI :: Magnitude');
                ismrm_imshow(angle(im),[],[1 size(im,3)],titles,'Synthetic MRI :: Phase');
                
                
        disp('<-- DIXON_COIL_COMBINE');
        
        
        %% SNR Estimation
        if(DEBUG == 2)
            disp('--> SNR_ROI');
                disp('    :: SNR Water Area');
                SNR.we1 = SNR_ROI('Water Echo 1', im(:,:,1));
                SNR.we3 = SNR_ROI('Water Echo 3', im(:,:,3), SNR.we1.noise.mask, SNR.we1.signal.mask);
                str_msg = fprintf('        -- Echo 1 : %.2f\n        -- Echo 3 : %.2f\n', SNR.we1.val, SNR.we3.val);
                disp('    :: SNR Fat Area');
                SNR.fe1 = SNR_ROI('Fat Echo 1', im(:,:,1), SNR.we1.noise.mask);
                SNR.fe3 = SNR_ROI('Fat Echo 3', im(:,:,3), SNR.we1.noise.mask, SNR.fe1.signal.mask);
                str_msg = fprintf('        -- Echo 1 : %.2f\n        -- Echo 3 : %.2f\n', SNR.fe1.val, SNR.fe3.val);
             disp('<-- SNR_ROI');
        end
         
        %% 3 points Dixon reconstruction
        if(DEBUG == 1 || DEBUG == 2)
            disp('--> DIXON_3P');
                disp('    :: Dixon processing');
                    [W, F] = Dixon_3P(im(:,:,1),im(:,:,2),im(:,:,3),TE);

                disp('    :: Water and Fat images');
                figure('Name','Synthetic MRI :: Dixon');
                subplot(121); imagesc(abs(W)); colormap gray, axis image; title('Water');
                subplot(122); imagesc(abs(F)); colormap gray, axis image; title('Fat');
            disp('<-- DIXON_3P');
        end
        
        %% Save data
         clear acc_factor ans data i id img L_normal_7 ne noise_color noise_white Rn_normal_7 s S0_F S0_W x y z
       % save('Synthetic.mat');
disp('<-- GENERATE_NOISY_MULTICOILS_FLASH_SIMULATION');