%%%%                        Script to get a csm simulation
%%%%
%%%% Protocole : 
%%%%        1. Generate a .h5 from /usr/local/bin/ismrmrd_generate_cartesian_shepp_logan
%%%%           with matrix dim, ncoils, noise desired (matrix dim should
%%%%           match encX, encY from Generate_Synthetic_FLASH_image)
%%%%        2. Extract the csm from the .h5 file
%%%%        3. Simulate FLASH image with ncoils with/without noise
%%%%        4. Go to ../../Matlab/FID2ccm and apply the DIXON

disp('--> GENERATE_NOISY_MULTICOILS_FLASH_SIMULATION');
    disp('    :: Dixon Fat/Water application');
        addpath('../../Matlab/');
        addpath('../../Matlab/CSM');

        %% Load session data
        [~,id]= system('whoami');
        str_user= id(1:end-1);
        %check_if_iam_using_the_ihuserver(str_user);
        %[ str_network_imagerie, str_network_perso ] = get_network_name( str_user );
        
        clear id

        %% Load the filename and data
        %filename = ['/home/', str_user, '/Dev/PostProcessing/Matlab/example.h5'];
        filename = 'example.h5';

        if exist(filename, 'file')
            hinfo = hdf5info(filename);   
            smaps = h5read(filename,  hinfo.GroupHierarchy.Groups(1).Datasets(2).Name);  
            smaps = complex(smaps.real, smaps.imag);
            ncoils = size(smaps,3);
        else
            error(['File ' filename ' does not exist.  Please generate it.'])
        end  
        
        Generate_Synthetic_FLASH_images;

        % Check if simulation is 2D
        if(size(synFLASH,3) == 1)
            synFLASH = permute(synFLASH, [1 2 4 3]);
        end        

        %% Use of ismrm_hansen_generate_data.m by Michael S. Hansen but softly modified
    disp('    :: Noise application & 7 channels simulation');
        addpath('../../ismrm_sunrise_matlab/');

        load('noise_covariances.mat');

        Rn_normal_7 = Rn_normal_8(1:ncoils,1:ncoils);
        Rn_broken_7 = Rn_broken_8(1:ncoils,1:ncoils);

        clear Rn_normal_8 Rn_broken_8

        L_normal_7 = chol(Rn_normal_7,'lower');
        L_broken_7 = chol(Rn_broken_7,'lower');

        %Some settings
        acc_factor = 1;

        for i=1:ne
            img = synFLASH(:,:,i);
            noise_level = 0.05*max(img(:));

            %Sample Cartesian Data
            noise_white = noise_level*complex(randn(size(img,1),size(img,2),size(smaps,3)),randn(size(img,1),size(img,2),size(smaps,3)));
            noise_color = reshape(permute(L_broken_7 * permute(reshape(noise_white, numel(noise_white)/ncoils,ncoils),[2 1]),[2 1]),size(noise_white));

            [data, sp] = ismrm_sample_data(img, smaps, acc_factor, 0);
            synData(:,:,:,i) = data + noise_color .* repmat(sp > 0,[1 1 ncoils]);
        end

        %% Coils combination - Image reconstruction
        disp('--> DIXON_COIL_COMBINE');
            disp('    :: Inati & Walsh reconstruction');

                % Number of echoes
                nechoes = 3;
                [im.inati, im.walsh] = Dixon_coil_combine(synData, nechoes, 'h5');
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
        
        %% Save data
        clear acc_factor ans data i id img L_normal_7 L_broken_7 ne noise_color noise_white Rn_broken_7 Rn_normal_7 s S0_F S0_W x y z
       % save('Synthetic.mat');
disp('<-- GENERATE_NOISY_MULTICOILS_FLASH_SIMULATION');