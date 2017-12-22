function [kspace_data, output_tmp] = hdf5_kspace_reader(path, files)

    [status,id]= system('whoami');

    str_user= id(1:end-1);

    check_if_iam_using_the_ihuserver(str_user);

    [ str_network_imagerie, str_network_perso ] = get_network_name( str_user );
    
    fprintf('\t Number of echoes : %d \n', size(files, 2));

    for ne=1:size(files, 2)

        number=num2str(files(ne));
        
        filename = ['/home/', str_user, path ,number,'.h5'];
        
        if exist(filename, 'file')
            dset = ismrmrd.Dataset(filename, 'dataset');
            fprintf('\t\t -- Echo %d : %s \n',ne, filename); 
        else
            error(['File ' filename ' does not exist.  Please generate it.'])
        end


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
            nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 ;
        catch
            nContrasts = 1;
        end

        D = dset.readAcquisition();



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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        contrast=1;
        slice=1;
        rep=1;

        % donne le nombre de ligne correspondant à ces parametres
        acqs_image = find(  (meas.head.idx.contrast==(contrast-1)) ...
            & (meas.head.idx.repetition==(rep-1)) ...
            & (meas.head.idx.slice==(slice-1))  );

        str_msg=sprintf('le nombre de lignes ACQ_IS_ ALL est  %d \n', size(acqs_image,2)); %disp(str_msg);

        mat.kspace.raw = zeros(enc_Nx, enc_Ny,  nCoils);

        for p = 1:length(acqs_image)

            ky = meas.head.idx.kspace_encode_step_1(acqs_image(p)) + 1;
            str_msg=sprintf('p %d  acqs_image(p)  %d  ky  %d', p, acqs_image(p), ky);  %disp(str_msg);

            echo = meas.head.idx.contrast(acqs_image(p)) + 1;

            mat.kspace.raw(:,ky,:) = meas.data{acqs_image(p)};

        end


        kspace_data(:,:,:,ne)= mat.kspace.raw;

    end

    output_tmp = ['/home/', str_user, path ,num2str(files(1)),'_',num2str(files(2)),'_',num2str(files(3))];