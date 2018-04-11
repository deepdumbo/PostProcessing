
%% Working with an existing ISMRMRD data set

% This is a simple example of how to reconstruct images from data
% acquired on a fully sampled cartesian grid
%
% Capabilities:
%   2D
%   use noise scans (pre-whitening)
%   remove oversampling in the readout direction
%   virtual coil using pca
%   coil reduction
%   magnitude and phase reconstruction
%   read data output from gadgetron
%
% Limitations:
%   only works with a single encoded space
%   fully sampled k-space (no partial fourier or undersampling)
%   one repetitions
%   doesn't handle averages, phases, segments and sets
%
%

% We first create a data set using the example program like this:
%   ismrmrd_generate_cartesian_shepp_logan -r 5 -C -o shepp-logan.h5
% This will produce the file shepp-logan.h5 containing an ISMRMRD
% dataset sampled evenly on the k-space grid -128:0.5:127.5 x -128:127
% (i.e. oversampled by a factor of 2 in the readout direction)
% with 8 coils, 5 repetitions and a noise level of 0.5
% with a noise calibration scan at the beginning


clear all

addpath('/home/valery/Reseau/Valery/MatlabUnix/ismrm_sunrise_matlab-master/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading an existing file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% filename_noise = '/home/valery/DICOM/20170124_Flash/FID/00112_NOISEFID01336_fl_grappa2_integre.h5';
filename = '/home/valery/DICOM/testdata.h5';
filename_noise = '/home/valery/DICOM/testdata_noise.h5';

if exist(filename_noise, 'file')
    dset_noise = ismrmrd.Dataset(filename_noise, 'dataset');
    hdr_noise = ismrmrd.xml.deserialize(dset_noise.readxml);
    D_noise = dset_noise.readAcquisition();
else
    error(['File ' filename_noise ' does not exist.  Please generate it.'])
end


if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

check_agreement_with_gadgetron='Y';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them


hdr = ismrmrd.xml.deserialize(dset.readxml);


%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

% Number of slices, coils, repetitions, contrasts etc.
% We have to wrap the following in a try/catch because a valid xml header may
% not have an entry for some of the parameters

try
    nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    nSlices = 1;
end

try
    nCoils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    nCoils = 1;
end

try
    nReps = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    nReps = 1;
end

try
    nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 + 1;
catch
    nContrasts = 1;
end


% TODO add the other possibilites

%% Read all the data
% Reading can be done one acquisition (or chunk) at a time,
% but this is much faster for data sets that fit into RAM.

D = dset.readAcquisition();



% TODO add a pre-whitening example
% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan > 1
    noise = D.select(1:firstScan-1);
else
    noise = [];
end

meas  = D.select(firstScan:D.getNumber);
clear D;

nbLines=size(meas.data,2);

%nb de lignes
str_msg=sprintf('le nombre de lignes TOTAL est %d \n', nbLines); disp(str_msg);
 % on va donc prendre juste l'image

 
contrast=1;
rep=1;
slice=1;

% donne le nombre de ligne qui correspondent à ces valeurs
acqs = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)));
 
 
K_image = zeros(enc_Nx, enc_Ny, enc_Nz,  nCoils);

size(K_image)

for p = 1:length(acqs)
    ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
    kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
    K_image(:,ky,kz,:) = meas.data{acqs(p)};
end



K1 = fftshift(ifft(fftshift(K_image,1),[],1),1);

K2 = fftshift(ifft(fftshift(K1,2),[],2),2);

im = fftshift(ifft(fftshift(K2,3),[],3),3);

center_z=size(im,3)/2;

figure(10)
for c = 1:nCoils
    subplot(2,4,c); imagesc(abs(im(:,:,center_z,c)));
end



addpath('/home/valery/Reseau/Valery/MatlabGit/Reconstruction/')

%% Take noise scans for pre-whitening example

% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isNoise = D_noise.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');

% toutes les données sont du bruit

firstScan = find(isNoise==1,1,'first');
if firstScan > 1
    noise = D_noise.select(1:firstScan-1);
else
    noise = [];
end

meas_noise  = D_noise.select(firstScan:D_noise.getNumber);
% clear D_noise;

%% pre-whitening (NoiseAdjustGadget.cpp NoiseAdjustGadget.h)

nbLines_noise=size(meas_noise.data,2);

%nb de lignes
str_msg=sprintf('le nombre TOTAL de lignes dans le fichier noise est %d \n', nbLines_noise); disp(str_msg);


%% calcul de la matrice de covariance
% pour commencer nous allons extraire les lignes correspondants au bruit
% ouvre l'ADC et acquiert le bruit electronique

% la première étape de traitement consiste à calculer la matrice de
% covariance du bruit
% cf fonction suivant du cours ismrm sunrise course 2013
% dmtx = ismrm_calculate_noise_decorrelation_mtx(noise(sp>0,:));

acq_noise_measurement =  find( (meas_noise.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT')) & ~(meas_noise.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA')) );
str_msg=sprintf('le nombre TOTAL de lignes dans noise  est %d \n', size(acq_noise_measurement,2)); disp(str_msg);

% nous allons utiliser 256 lignes

% allocation de la matrice de covariance du bruit
% la taille correspondant aux nombres d'antennes



% pour chaque ligne, en tenant compte de l'oversampling, nous disposons
% d'une matrice à deux dimensions [enc_Nx, nCoils ]
% nous effections une permutation
% nous calculons le produit du bruit et de sa transposée
% nou sommons les données et comptons le nombre de points utilisés

[ matlab.covariance_matrix.triangular ,  matlab.covariance_matrix.raw_normalise, matlab.covariance_matrix.raw] = extract_covariance_matrix( meas_noise, acq_noise_measurement , nCoils );

[ data_pre_whitening ,  matlab.covariance_matrix.triangular_bw  , data_before_pre_whitening] = apply_pre_whitening_noise( meas, meas_noise, hdr, matlab.covariance_matrix.triangular, 0 );

data_pre_whitening_ok=reshape(data_pre_whitening, [128, 128, 128, 4 ]);


Image_pre_whitening = zeros(enc_Nx, enc_Ny, enc_Nz,  nCoils);



figure(12) 
subplot(121); imagesc(abs(matlab.covariance_matrix.raw)) ; title ('covariance matrix');
subplot(122); imagesc(abs(matlab.covariance_matrix.triangular)) ; title ('covariance matrix after triangulation');

K1 = fftshift(ifft(fftshift(data_pre_whitening_ok,1),[],1),1);

K2 = fftshift(ifft(fftshift(K1,2),[],2),2);

Image_pre_whitening = fftshift(ifft(fftshift(K2,3),[],3),3);


center_z=ceil(size(im,3)/3);

figure(11)
for c = 1:nCoils
     subplot(2,4,c); imagesc(abs(im(:,:,center_z,c))); title (['raw data, coil: ' num2str(c)]);
    subplot(2,4,c+nCoils); imagesc(abs(Image_pre_whitening(:,:,center_z,c))); title (['noise prewigtening, coil: ' num2str(c)]);
end



acqs_center = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & (meas.head.idx.kspace_encode_step_2==63));



