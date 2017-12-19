%% Reconstruction of 3D radial vibe data
clear all; close all; clc; clear classes;

addpath(genpath('./utils'));
addpath(genpath('../../gpuNUFFT'));
addpath('/home/valery/Reseau/Valery/MatlabUnix/ismrm_sunrise_matlab-master/');

number_of_output_respiratory_phases=3;
number_of_output_cardiac_phases=20;

% stop à 2 method , 2 resp , 2 cardiac

for method=1:2
    
    for j=1:1:number_of_output_respiratory_phases
        for i=1:number_of_output_cardiac_phases
            
            clear trajectory
            clear kspace
            
            str_i=sprintf('%d',i );disp(str_i);
            str_j=sprintf('%d',j );disp(str_j);
            
            filename_traj=['/home/valery/Tempo/' , 'trajectory_subsets_', str_j , '_', str_i , '.mat'];
            filename_kspace=['/home/valery/Tempo/' , 'kspace_subsets_', str_j , '_', str_i , '.mat'];
            
            A=load(filename_kspace);
            B=load(filename_traj);
            
            kspace_subsets=A.kspace;
            trajectory_subsets=B.trajectory;
            
            size(kspace_subsets);
            size(trajectory_subsets);             
            
            disp('starting density compensation');
            info.verbose = 0;
            info.numIter = 1;
            info.effMtx  = size(A.kspace,1);
            info.osf=1;
            DCF = sdc3_MAT(reshape(trajectory_subsets,[],3)',info.numIter,info.effMtx,info.verbose,info.osf);
            
            
            
            %% Data parameters
            info.N=size(kspace_subsets,1);
            info.nSl=info.N;
            nCh=size(kspace_subsets,3);
            disp_slice=info.N/2;
            useGPU = true;
            useMultiCoil = 1;
            
            
            %% Load data
            % load ./data/rawdata_phantom_regridding.mat;
            % [nPE,nFE,nCh]=size(rawdata);
            % rawdata = reshape(rawdata,[nPE*nFE,nCh]);
            
            %% Regridding operator GPU without coil sensitivities for now
            disp('Generate NUFFT Operator without coil sensitivities');
            osf = 2; wg = 3; sw = 8;
            imwidth = info.N;
            
            % FT = gpuNUFFT(k',col(w(:,:,1)),osf,wg,sw,[N,N,nSl],[],true);
            info.wg =2;
            info.sw = 4;
            mattest=reshape(trajectory_subsets,[],3)';
            FT = gpuNUFFT(mattest,DCF,info.osf,info.wg,info.sw,[info.N,info.N,info.nSl],[],true);
            
            for ii=1:nCh
                img_sens(:,:,:,ii) = FT'*reshape(kspace_subsets(:,:,ii),[],1);
            end
            
            %% Estimate sensitivities
            disp('Estimate coil sensitivities.');                      
            
            senseEst=zeros(size(img_sens));
            
            if (method==1)
                
                % Use this instead for more reasonable sensitivitities, but takes some time
                for ii=1:info.nSl
                    if(mod(ii,10)==1)
                    disp(['Slice ', num2str(ii), '/', num2str(info.nSl)]);
                    end
                    [~,senseEst(:,:,ii,:)]=adapt_array_2d(squeeze(img_sens(:,:,ii,:)));
                end
                
            elseif (method==2)
                
                % use this instead
                for ii=1:info.nSl
                    if(mod(ii,10)==1)
                    disp(['Slice ', num2str(ii), '/', num2str(info.nSl)]);
                    end
                    [ senseEst(:,:,ii,:) ] = coil_map_study_2d_Inati(squeeze(img_sens(:,:,ii,:)), 5, 3 );
                end                
                
            else  
                                
                % Terribly crude, but fast
                 img_sens_sos = sqrt(sum(abs(img_sens).^2,4));
                 senseEst = img_sens./repmat(img_sens_sos,[1,1,1,nCh]);
            end
            
            
            %% Redefine regridding operator GPU including coil sensitivities
            disp('Generate NUFFT Operator with coil sensitivities');
            % FT = gpuNUFFT(k',col(w(:,:,1)),osf,wg,sw,[N,N,nSl],senseEst,true);
            
            FT = gpuNUFFT(mattest,DCF,info.osf,info.wg,info.sw,[info.N,info.N,info.nSl],senseEst,true);
            
            %% Forward and adjoint transform
            tic
            img_comb = FT'*reshape(kspace_subsets,[],nCh);
            timeFTH = toc;
            disp(['Time adjoint: ', num2str(timeFTH), ' s']);
            % figure,imshow(abs(img_comb(:,:,disp_slice)),[]); title('Regridding');
            % figure,kshow(abs(fft2c(img_comb(:,:,disp_slice)))); title('Regridding k-space');
            
            tic
            test = FT*img_comb;
            timeFT = toc;
            disp(['Time forward: ', num2str(timeFT), ' s']);
            
            
            %% Reconstruction parameters
            maxitCG = 20;
            alpha = 1e-5;
            tol = 1e-6;
            display = 1;
            
            %% CGSENSE Reconstruction
            mask = 1;
            tic
            img_cgsense = cg_sense_3d(reshape(kspace_subsets,[],nCh),FT,senseEst,mask,alpha,tol,maxitCG,display,disp_slice,useMultiCoil);
            timeCG = toc;
            disp(['Time CG SENSE: ', num2str(timeCG), ' s']);
                        
            if (method==1)
                str_method='adapt';
            elseif (method==2)
                str_method='inati';
            else
                str_method='';
            end
            
            filename_img_cgsense=['/home/valery/Tempo/' , 'img_cgsense_', str_j , '_', str_i , '_',str_method , '.mat'];
            filename_img_comb=['/home/valery/Tempo/' , 'img_comb_', str_j , '_', str_i , '_',str_method ,  '.mat'];
                        
            save(filename_img_comb,'img_comb');
            save(filename_img_cgsense,'img_cgsense');            
            
        end
    end    
end
    


    % disp_slice=100
    %
    % %% Display
    % figure;
    % subplot(1,2,1); imshow(abs(img_comb(:,:,disp_slice)),[]); title('Regridding');
    % subplot(1,2,2); imshow(abs(img_cgsense(:,:,disp_slice)),[]); title('CGSENSE');    
    % load('/home/valery/Tempo/img_comb_adapt.mat','img_comb');
    % load('/home/valery/Tempo/img_cgsense_adapt.mat','img_cgsense');