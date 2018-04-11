
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

function [mat] = demo_flash_grappa_generic_kylian(h5_path)
   % clear all

    addpath('../../Generic_functions/');

    [ str_user ] = get_PC_name();


    addpath(['/home/',str_user,'/Dev/PostProcessing/ismrm_sunrise_matlab/']);
    addpath('../');
    addpath('Utility/');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading an existing file %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filename = [ h5_path,'.h5'];
    if exist(filename, 'file')
        dset = ismrmrd.Dataset(filename, 'dataset');
    else
        error(['File ' filename ' does not exist.  Please generate it.'])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Read some fields from the XML header %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We need to check if optional fields exists before trying to read them

    % hdr_noise = ismrmrd.xml.deserialize(dset_noise.readxml);
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

     [ number_of_slices, number_of_channels  , number_of_repetitions, number_of_contrasts, number_of_phase, number_of_average , number_of_segments] = get_number_of( hdr );

    if (number_of_contrasts>1)

        disp(' attention si multiecho inati ne fonctionne pas car tu appliques à chauqe fois des coefficients differents');

    end


    % TODO add the other possibilites

    %% Read all the data
    % Reading can be done one acquisition (or chunk) at a time,
    % but this is much faster for data sets that fit into RAM.
    D = dset.readAcquisition();

    %% Take noise scans for pre-whitening example

    %% quelques vérifications

    if (hdr.encoding.encodedSpace.matrixSize.x~=hdr.encoding.reconSpace.matrixSize.x)

        disp(' presence d oversampling \n');

    end

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

    size(meas.data)

    for p = 1:length(meas.data)

        data_before_pre_whitening(:,p,:)=double(meas.data{p});

    end


    %% do the reco % do a loop for slice and contrast

    contrast=1;
    slice=1;
    rep=1;
    set=1;

    % donne le nombre de ligne correspondant à ces parametres
    acqs_image_all = find(  (meas.head.idx.contrast==(contrast-1)) ...
        & (meas.head.idx.repetition==(rep-1)) ...
        & (meas.head.idx.slice==(slice-1))...
        & (meas.head.idx.set==(set-1))   );

%     str_msg=sprintf('le nombre de lignes ACQ_IS_ ALL est  %d \n', size(acqs_image_all,2)); disp(str_msg);

    % on passe les donnees en kspace

    mat.kspace.raw = zeros(enc_Nx, enc_Ny,  number_of_channels);
    mat.data = zeros(enc_Nx, enc_Ny);


    for p = 1:length(acqs_image_all)

        ky = meas.head.idx.kspace_encode_step_1(acqs_image_all(p)) + 1;
%         str_msg=sprintf('p %d  acqs_image_all(p)  %d  ky  %d', p, acqs_image_all(p), ky);  disp(str_msg);

        mat.kspace.raw(:,ky,:) = meas.data{acqs_image_all(p)};
        mat.data(:,ky)=1;

    end




%     figure()
%     imagesc(mat.data)


    % donne le nombre de ligne correspondant à ces parametres
    acqs_image_only = find(  (meas.head.idx.contrast==(contrast-1)) ...
        & (meas.head.idx.repetition==(rep-1)) ...
        & (meas.head.idx.slice==(slice-1))...
        & (meas.head.idx.set==(set-1))...    
        & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

%     str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE est  %d \n', size(acqs_image_only,2)); disp(str_msg);

    % donne le nombre de ligne correspondant à ces parametres
    acqs_parallel_calibration = find(  (meas.head.idx.contrast==(contrast-1)) ...
        & (meas.head.idx.repetition==(rep-1)) ...
        & (meas.head.idx.slice==(slice-1))...
        & (meas.head.idx.set==(set-1)) ...
        & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

%     str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION est  %d \n', size(acqs_parallel_calibration,2)); disp(str_msg);

    % donne le nombre de ligne correspondant à ces parametres
    acqs_parallel_calibration_and_imaging = find(  (meas.head.idx.contrast==(contrast-1)) ...
        & (meas.head.idx.repetition==(rep-1)) ...
        & (meas.head.idx.slice==(slice-1))...
        & (meas.head.idx.set==(set-1)) ...
        & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING'))  );

%     str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING est  %d \n', size(acqs_parallel_calibration,2)); disp(str_msg);

    %% gadget remove oversampling (RemoveROOversamplingGadget.cpp)

    % ce gadget effectue une transformée de fourrier inverse dans la direction x
    % on crop
    % on fait une transformée de fourrier dans la direction x


    %% partie PCA and coil compression (PCA coils gadgets)



    %% nous allons maintenant attaquer les gadgets GrappaGadget GrappaUnimixingGadget

    % le kspace est strictement identique entre maltab et gadgetron ?

    sp = zeros(enc_Nx, enc_Ny);

    for p = 1:length(acqs_image_only)
        ky = meas.head.idx.kspace_encode_step_1(acqs_image_only(p)) + 1;
        sp(:,ky) = 1;
    end

    for p = 1:length(acqs_parallel_calibration)
        ky = meas.head.idx.kspace_encode_step_1(acqs_parallel_calibration(p)) + 1;
        sp(:,ky) = sp(:,ky)+ 2;
    end

    for p = 1:length(acqs_parallel_calibration_and_imaging)
        ky = meas.head.idx.kspace_encode_step_1(acqs_parallel_calibration_and_imaging(p)) + 1;
        sp(:,ky) = sp(:,ky)+ 2;
    end


%     imagesc(sp)




    [mat.image.grappa.walsh] = ismrm_cartesian_GRAPPA(mat.kspace.raw(:,:,:),sp,2);

%     [mat.image.grappa.walsh_kellman] = ismrm_cartesian_GRAPPA_val(mat.kspace.raw(:,:,:),sp,2, [], [] , 1);

%     [mat.image.grappa.walsh_python] = ismrm_cartesian_GRAPPA_val(mat.kspace.raw(:,:,:),sp,2, [], [] , 2);

%     [mat.image.grappa.inati_matlab] = ismrm_cartesian_GRAPPA_val(mat.kspace.raw(:,:,:),sp,2, [], [] , 3);



%     figure(4)
%     subplot(2,4,1);  imagesc(abs(mat.image.grappa.walsh_kellman)); 
%     subplot(2,4,2);  imagesc(angle(mat.image.grappa.walsh_kellman)); colormap(gray); title('walsh 1');
%     subplot(2,4,3);  imagesc(abs(mat.image.grappa.walsh_python));
%     subplot(2,4,4);  imagesc(angle(mat.image.grappa.walsh_python)); colormap(gray);  title('walsh 2 ');
%     subplot(2,4,5);  imagesc(abs(mat.image.grappa.inati_matlab));
%     subplot(2,4,6);  imagesc(angle(mat.image.grappa.inati_matlab)); colormap(gray);  title('inati');
%     subplot(2,4,7);  imagesc(abs(mat.image.grappa.walsh));
%     subplot(2,4,8);  imagesc(angle(mat.image.grappa.walsh)); colormap(gray);  title('walsh');
end