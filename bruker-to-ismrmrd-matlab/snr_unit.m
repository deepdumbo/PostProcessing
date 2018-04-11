clear all
[ str_user ] = get_PC_name();

% [ str_network_imagerie, str_network_perso ] = get_network_name( str_user );

filename=['/home/',str_user,'/Dev/Data/out_snr_unit.h5'];
hinfo = hdf5info(filename);

% hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name

res_1=single(h5read(filename, hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name));    
res_2=single(h5read(filename, hinfo.GroupHierarchy.Groups(1).Groups(2).Datasets(2).Name));









filename_noise='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-24_Flash/FID/00110_NOISEFID01334_fl_simple.h5'
filename_data='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-24_Flash/FID/00110_FID01334_fl_simple.h5'

% filename_noise = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-24_Flash/FID/00112_NOISEFID01336_fl_grappa2_integre.h5';
% filename_data = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-24_Flash/FID/00112_FID01336_fl_grappa2_integre.h5';


if exist(filename_noise, 'file')
    dset_noise = ismrmrd.Dataset(filename_noise, 'dataset');
else
    error(['File ' filename_noise ' does not exist.  Please generate it.'])
end


if exist(filename_data, 'file')
    dset = ismrmrd.Dataset(filename_data, 'dataset');
else
    error(['File ' filename_data ' does not exist.  Please generate it.'])
end


hdr_noise = ismrmrd.xml.deserialize(dset_noise.readxml);
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



% TODO add the other possibilites

%% Read all the data
% Reading can be done one acquisition (or chunk) at a time,
% but this is much faster for data sets that fit into RAM.
D_noise = dset_noise.readAcquisition();
D = dset.readAcquisition();




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
clear D_noise;




%% Ignore noise scans
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
% clear D;




acq_noise_measurement =  find( (meas_noise.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT'))  ...
                                & ~(meas_noise.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA')) );

str_msg=sprintf('le nombre TOTAL de lignes dans l acquisition de bruit  est %d \n', size(acq_noise_measurement,2)); disp(str_msg);
str_msg=sprintf('le nombre TOTAL de lignes dans le fichier image est %d \n', size(meas.data,2)); disp(str_msg);

[ triangular_covariance_matrix ] = extract_covariance_matrix( meas_noise, acq_noise_measurement , nCoils );

oversampling=0;
[ data_pre_whitening ,  triangular_covariance_matrix_bw  , data_before_pre_whitening] = apply_pre_whitening( meas, meas_noise, hdr, triangular_covariance_matrix , oversampling );

contrast=1;
slice=1;
rep=1;
set=1;

% donne le nombre de ligne correspondant à ces parametres
acqs_image = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & (meas.head.idx.set==(set-1))   );

str_msg=sprintf('le nombre de lignes ACQ_IS_ ALL est  %d \n', size(acqs_image,2)); disp(str_msg);


for p = 1:length(acqs_image)
        
        ky = meas.head.idx.kspace_encode_step_1(acqs_image(p)) + 1;        
        str_msg=sprintf('p %d  acqs_image_all(p)  %d  ky  %d', p, acqs_image(p), ky);  disp(str_msg);        
        str_e1=sprintf('%d',ky-1);
        
        kspace_oversampling(:,ky,:)=data_pre_whitening(:,acqs_image(p),:);                
        
end


image_oversampling=ifft_2D(kspace_oversampling);

image=image_oversampling(65:64+128,:,:);

csm_inati = coil_map_study_2d_Inati( image, 5, 3 );

csm_inati = ismrm_normalize_shading_to_sos(csm_inati);

ccm_inati = ismrm_compute_ccm(csm_inati);

res_3 = abs(sum(image .* ccm_inati, 3));



figure(1)
subplot(131);imagesc(res_1)
subplot(132);imagesc(res_2)
subplot(133);imagesc(res_3)


max(res_1(:))/max(res_3(:))
min(res_1(:))/min(res_3(:))



figure(1)
subplot(131);imagesc(res_1)
subplot(132);imagesc(res_3*161)
subplot(133);imagesc(res_1-res_3*161)