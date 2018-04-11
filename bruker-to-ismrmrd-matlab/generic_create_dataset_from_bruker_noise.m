%% Read Bruker File

clear all
close all


%% TODO voici la liste des problemes
% completer les header
% * si 3D encoding limit il ne faut peut -être pas envoyer la dernier ligne
% * si acceleration factor y different de 2
% * si acceleration factor z
% * si average > 1
% * si echo > 1
% * si stacks > 1
% * si partial fourier y

%% voici les reco qui fonctionne 
% * E5 2D normal  (128*128)
% * E6 2D GRAPPA Y 1 Slice  (76*128)
% * E7 2D 3 slices  (128*128*3)
% * E8 2D GRAPPA Y 3 Slices (76*128*3)
% * E9 = E8 2 repetitions
% * E10 =E6 2 repetitions
% * E11 3D GRAPPA Y
% * E12 3D normal
% * E72 2D complet (matrix size 128*128)  
% * E73 2D complet (matrix size 114*114) 
% * E74 2D complet (matrix size 99*100) 
% * E75 2D partial fourier (matrix size 128*128) 
% * E76 2D partial fourier (matrix size 114*114) 
% * E77 2D partial fourier (matrix size 99*100) 

%% voici les reco qui ne fonctionne pas
% *mp2rage

% acquisition_path='/home/valery/DICOM/badres/4/'
% acquisition_path='/home/valery/DICOM/Bruker_7617/11/';
% acquisition_path='/home/valery/DICOM/tests13062017/80/';
% acquisition_path='/home/valery/DICOM/tests13062017/74/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Sheep/Control/2014_04_10_Heart_1/77/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Kadence/Infarct/2015_12_16_Heart_3/31/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Sheep/Control/2015_12_13_Heart_4/31/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Kadence/Control/2016_05_03_Heart_2/17/';
% acquisition_path='/home/valery/DICOM/Liryc/noise/20170922_094726_AT_Liryc_Noise01_1_1/3/';
% acquisition_path='/home/valery/DICOM/Liryc/noise/20170922_094726_AT_Liryc_Noise01_1_1/6/';
output_filename = '/home/valery/DICOM/testdata_bruker_noise.h5';
acquisition_path='/home/valery/DICOM/Noise/8/';


%% reading bruker acqp, method and fid files. 

ex = read_bru_experiment(acquisition_path);
% ex = read_bru_experiment_without_fid(acquisition_path);
% 
% % read the raw data
% if exist([acquisition_path '/fid'])
%     fpre = fopen([acquisition_path '/fid'],'r');
%     fid = fread(fpre,'int32','l');
%     fid = fid(1:2:end) - 1i * fid(2:2:end);
%     fclose(fpre);
% end
% 
% ex.fid=single(fid);
% clear fid


% load('/home/valery/DICOM/Noise/8/LignesAcq_noheader.txt');
tempo=readParams_Bruker2( 'dirPath',   acquisition_path);

A=tempo.LignesData;

size(A)

B=reshape(A, [], 2);

size(B)

C=complex(B(:,1), B(:,2));

noise_from_acqp=reshape(C, 128, 7, []);

size(noise_from_acqp)


[ nX, nY, nZ ] = get_dimensions( ex );

[ nX_acq, nY_acq, nZ_acq ] = get_dimensions_acq( ex );

[ readout, E1, E2 ] = get_encoding_size( ex , nZ );

[header] = fill_the_flexible_xml_header_noise(ex);

%% reshape the fid to match the fixed data structure and remove the zero

% [ data_for_acqp ] = remove_zero_from_fid( ex );

number_of_channels=size(noise_from_acqp,2);
number_of_lines=size(noise_from_acqp,3);

%% gestion de la diffusion

% size(data_for_acqp)
% 
% data_tempo=reshape(data_for_acqp,[readout,number_of_channels , 10 ,E1 , E2 ]);
% 
% ex.method.PVM_EncSteps1(end)-ex.method.PVM_EncSteps1(1)
% 
% data_ok=zeros(size(data_tempo,1),size(data_tempo,2),size(data_tempo,3),nY,size(data_tempo,5));
% 
% data_ok(:,:,:,ex.method.PVM_EncSteps1(1)+round(nY/2):ex.method.PVM_EncSteps1(end)+round(nY/2),:)=data_tempo;
% 
% figure()
% imagesc(abs(squeeze(data_ok(:,1,1,:,75))))
% 
% data_reco=permute(data_ok, [1, 4, 5, 2,3]);
% 
% size(data_reco)
% 
% img_reco=ifft_3D(data_reco);
% 
% figure()
% for s=20:100
%     for d=1:10
% subplot(2,5,d); imagesc(abs(squeeze(img_reco(:,:,s,1,d))), [0 1]); colormap(gray);
% end
% pause(0.2)
% end
% 
% ex.acqp.ACQ_time_points
% ex.method.PVM_DwDir
% ex.method.PVM_DwNDiffExp
%% fill the fixed data structure

[ idx, flag ] = fill_the_idx_noise( header , ex);


%% parallel imaging option

[ Ysamp_regular , Ysamp_ACS , acceleration_factor_y] = check_options_for_parallel_imaging( ex );


%% Generating a simple ISMRMRD data set

% This is an example of how to construct a datset from synthetic data
% simulating a fully sampled acquisition on a cartesian grid.
% data from 4 coils from a single slice object that looks like a square

% File Name

delete(output_filename)
% Create an empty ismrmrd dataset
if exist(output_filename,'file')
    error(['File ' output_filename ' already exists.  Please remove first'])
end
dset = ismrmrd.Dataset(output_filename);

% It is very slow to append one acquisition at a time, so we're going
% to append a block of acquisitions at a time.
% In this case, we'll do it one repetition at a time to show off this
% feature.  Each block has nY aquisitions
acqblock = ismrmrd.Acquisition(number_of_lines);

% Set the header elements that don't change
acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = readout;
acqblock.head.center_sample(:) = floor(readout/2);
acqblock.head.active_channels(:) = number_of_channels;
acqblock.head.read_dir  = repmat([1 0 0]',[1 number_of_lines]);
acqblock.head.phase_dir = repmat([0 1 0]',[1 number_of_lines]);
acqblock.head.slice_dir = repmat([0 0 1]',[1 number_of_lines]);

% Loop over the acquisitions, set the header, set the data and append

for acqno = 1:number_of_lines
    
    % Set the header elements that change from acquisition to the next
    % c-style counting
    acqblock.head.scan_counter(acqno) =  acqno-1;
    acqblock.head.idx.kspace_encode_step_1(acqno) = idx.kspace_encode_step_1(acqno);
    acqblock.head.idx.kspace_encode_step_2(acqno) = idx.kspace_encode_step_2(acqno);
    acqblock.head.idx.repetition(acqno) = idx.repetition(acqno);
    acqblock.head.idx.slice(acqno) = idx.slice(acqno);
    acqblock.head.idx.set(acqno) = idx.set(acqno);
    acqblock.head.idx.contrast(acqno) = idx.contrast(acqno);
    
    % Set the flags
    acqblock.head.flagClearAll(acqno);
    acqblock.head.sample_time_us(acqno)=1/ex.acqp.SW_h*1e6;
    
    %     str_msg=sprintf('count %d Ysamp %d        ', acqno,   idx.kspace_encode_step_1(acqno)    ); disp( str_msg);
    
%     if (acceleration_factor_y>1)
%         
%         if ismember(idx.kspace_encode_step_1(acqno),Ysamp_ACS-1)
%             %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne
%             
%             if ismember(idx.kspace_encode_step_1(acqno),Ysamp_regular)
%                 % both calibration and part of the undersampled pattern
%                 acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
%                 %                 disp('both calibration and part of the undersampled pattern');
%             else
%                 % in ACS block but not part of the regular undersampling
%                 acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
%                 %                 disp('in ACS block but not part of the regular undersampling');
%             end
%         end
%     end
    
     acqblock.head.flagSet('ACQ_IS_NOISE_MEASUREMENT', acqno)
    
    
    if (flag.first_in_encoding_step1(acqno)== 1)
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
    end
    
    if (flag.last_in_encoding_step1(acqno)== 1)
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
    end
    
    %         acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
    %         acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
    %     elseif acqno==size(K,2)
    %         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
    %         acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
    %         acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
    %     end    
    
    % fill the data
    acqblock.data{acqno} = squeeze(noise_from_acqp(:,:,acqno));
    
end

% Append the acquisition block
dset.appendAcquisition(acqblock);



%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%

header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = acceleration_factor_y ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
header.encoding.parallelImaging.calibrationMode = 'embedded' ;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();

fprintf('-------------------------------------------------------------- \n');
disp('fin de la conversion');
fprintf('-------------------------------------------------------------- \n');














number_of_channels=size(noise_from_acqp,2);

dmtxTotal=zeros(number_of_channels,number_of_channels);

Mtotal=0;

for p=1:size(noise_from_acqp,3);
    
    noise_samples=double(noise_from_acqp(:,:,p));
%     str_msg=sprintf('ligne %d  ,  matrix size [ %d , %d ] \n', p , size(noise_samples,1) , size(noise_samples,2));disp(str_msg);
    noise = reshape(noise_samples,numel(noise_samples)/size(noise_samples,length(size(noise_samples))), size(noise_samples,length(size(noise_samples))));
%     noise = permute(noise,[2 1]);
    
    M = size(noise,1);
    
    % Performs generalized matrix multiplication.
    dmtx = (noise'*noise);
       
    
    % nous comptons le nombre de points utilisés pour la normalisation
    % future
    Mtotal=Mtotal+M;
    dmtxTotal=dmtxTotal+dmtx;
    
end

str_msg=sprintf('le nombre de points utilisés est %d  ', Mtotal);disp(str_msg);

% normalisation
dmtxTotalNormalize = (1/(Mtotal-1))*dmtxTotal;

% Triangulation of matrices

% le facteur est important sinon cela ne correspond pas avec le Gadgetron
% TODO trouver une meilleur raison

scale_factor=1;

dmtxTotalNormalize_temp = tril(inv(chol(dmtxTotalNormalize,'lower')));
dmtxTotalNormalize_triangulation = dmtxTotalNormalize_temp*sqrt(2)*sqrt(scale_factor);




figure(12) 
subplot(121); imagesc(abs(dmtxTotalNormalize)) ; title ('covariance matrix');
subplot(122); imagesc(abs(dmtxTotalNormalize_triangulation)) ; title ('covariance matrix after triangulation');



disp('sample time');
disp(1/ex.acqp.SW_h*1e6)



% 
% kspace=permute(data_for_acqp, [1 3 2]);
% image_test=ifft_2D(kspace);
% magnitude=sqrt(sum(abs(image_test).^2,3));
% 
% figure(1)
% subplot(131); imagesc(magnitude(:,:)); colormap(gray);
% 
% 
% 
% %%  pre-whitening
% in=kspace;
% orig_size = size(in);
% nelements = prod(orig_size)/number_of_channels;
% 
% in = permute(reshape(in,nelements,number_of_channels),[2 1]);
% out = reshape(permute(dmtxTotalNormalize_triangulation*in,[2 1]),orig_size);
% 
% data_pre_whitening=out;
% 
% image_test_pre_whitening=ifft_2D(data_pre_whitening);
% magnitude_pre_whitening=sqrt(sum(abs(image_test_pre_whitening).^2,3));
% 
% subplot(132); imagesc(magnitude_pre_whitening(:,:)); colormap(gray);
% 
% subplot(133); imagesc(magnitude_pre_whitening(:,:)-magnitude(:,:) ); colormap(gray);

