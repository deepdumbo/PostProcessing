clear all


[status,id]= system('whoami');

str_user= id(1:end-1);

check_if_iam_using_the_ihuserver(str_user);

[ str_network_imagerie, str_network_perso ] = get_network_name( str_user );


% filename='/home/valery/generic_spirit_random_cine.h5'
filename=['/home/',str_user,'/DICOM/Echo/data_te8.h5'];

% filename='/home/valery/testdata.h5';
% filename='/home/valery/Reseau/Imagerie/For_Valery/testdata.h5'

if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

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
    number_of_slices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    number_of_slices = 1;
end

try
    number_of_coils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    number_of_coils = 1;
end

try
    number_of_repetitions = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    number_of_repetitions = 1;
end

try
    number_of_contrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 ;
catch
    number_of_contrasts = 1;
end

try
    number_of_sets = hdr.encoding.encodingLimits.set.maximum + 1 ;
catch
    number_of_sets = 1;
end

D = dset.readAcquisition();

% Note: can select a single acquisition or header from the block, e.g.
%  acq = D.select(1);
%  hdr = D.head.select(1);
% or you can work with them all together





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
clear D;



contrast=1;
rep=1;
set=1;


acqs_image_only = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.set==(set-1))...    
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE est  %d \n', size(acqs_image_only,2)); disp(str_msg);


contrast=1;
acqs_image_k1_k2_c1 = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.set==(set-1))... 
    & (meas.head.idx.kspace_encode_step_2==(enc_Nz/2))...
    & (meas.head.idx.kspace_encode_step_1==(enc_Ny/2))...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE est  %d \n', size(acqs_image_k1_k2_c1,2)); disp(str_msg);



contrast=2;
acqs_image_k1_k2_c2 = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.set==(set-1))... 
    & (meas.head.idx.kspace_encode_step_2==(enc_Nz/2))...
    & (meas.head.idx.kspace_encode_step_1==(enc_Ny/2))...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE est  %d \n', size(acqs_image_k1_k2_c2,2)); disp(str_msg);

clear tempo
tempo=meas.data(acqs_image_k1_k2_c1);
fid_c1=tempo{1,1};

clear tempo
tempo=meas.data(acqs_image_k1_k2_c2);
fid_c2=tempo{1,1};


figure(1)
plot(abs(fid_c1(:,1)));
hold on;
plot(abs(fid_c2(:,1)));



contrast=1;
acqs_image_c1 = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.set==(set-1))...   
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE_contrast1 est  %d \n', size(acqs_image_c1,2)); disp(str_msg);



 for p = 1:length(acqs_image_c1)
                ky = meas.head.idx.kspace_encode_step_1(acqs_image_c1(p)) + 1;
                kz = meas.head.idx.kspace_encode_step_2(acqs_image_c1(p)) + 1;
                mat.kspace.siemens(:,ky,kz,:)=meas.data{acqs_image_c1(p)};
               
 end
 
 
contrast=2;
acqs_image_c2 = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.set==(set-1))...   
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE_contrast1 est  %d \n', size(acqs_image_c2,2)); disp(str_msg);



 for p = 1:length(acqs_image_c2)
                ky = meas.head.idx.kspace_encode_step_1(acqs_image_c2(p)) + 1;
                kz = meas.head.idx.kspace_encode_step_2(acqs_image_c2(p)) + 1;
                mat.kspace.siemens2(:,ky,kz,:)=meas.data{acqs_image_c2(p)};
               
 end

 
mat.image.siemens = ifft_3D( mat.kspace.siemens);
mat.image.siemens2 = ifft_3D( mat.kspace.siemens2);


close(figure(1))
figure(1)

for s=1:size(mat.image.siemens,3)
subplot(231); imagesc(abs(mat.image.siemens(:,:,s,1)));  colormap(gray);
subplot(232); imagesc(abs(mat.image.siemens2(:,:,s,1)));  colormap(gray);
subplot(233); imagesc(abs(mat.image.siemens(:,:,s,1))-abs(mat.image.siemens2(:,:,s,1)));  colormap(gray);
subplot(234); imagesc(angle(mat.image.siemens(:,:,s,1)));  colormap(gray);
subplot(235); imagesc(angle(mat.image.siemens2(:,:,s,1)));  colormap(gray);
subplot(236); imagesc(angle(mat.image.siemens(:,:,s,1))-angle(mat.image.siemens2(:,:,s,1)));  colormap(gray);
pause(0.05)
end
