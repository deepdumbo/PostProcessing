%% Read Bruker File

clear all
close all

[status,id]= system('whoami');

% str_user= id(1:end-1);
str_user = get_PC_name();

check_if_iam_using_the_ihuserver(str_user);

[ str_network_imagerie, str_network_perso ] = get_network_name( str_user );

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

% acquisition_path= ['/home/', str_user, '/DICOM/Bruker_7617/11'];
% acquisition_path='/home/valery/DICOM/tests13062017/72/'
% acquisition_path='/home/valery/DICOM/badres/4/'
% acquisition_path='/home/valery/DICOM/Bruker_7617/11/';
% acquisition_path='/home/valery/DICOM/tests13062017/80/';
% acquisition_path='/home/valery/DICOM/tests13062017/74/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Sheep/Control/2014_04_10_Heart_1/77/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Kadence/Infarct/2015_12_16_Heart_3/31/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Sheep/Control/2015_12_13_Heart_4/31/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Kadence/Control/2016_05_03_Heart_2/17/';
% acquisition_path='/home/valery/DICOM/Liryc/noise/20170922_094726_AT_Liryc_Noise01_1_1/3/';
% acquisition_path='/home/valery/DICOM/Liryc/noise/20170922_094726_AT_Liryc_Noise01_1_1/3/';

% acquisition_path='/home/valery/Reseau/Imagerie/For_Kylian/Flash_Kylian/9';
% acquisition_path='/home/valery/DICOM/Noise/8/';

% acquisition_path='/home/valery/Reseau/Imagerie/For_Kylian/Test_comp_DixonVARPRO/29/'

numb = [9:11];

for i=numb(1):2:numb(end-2)
    acquisition_path_1 = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Verification/Ex_Vivo/2D/20180412_Unmasc/', num2str(i)];
    acquisition_path_2 = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Verification/Ex_Vivo/2D/20180412_Unmasc/', num2str(i+1)];
    acquisition_path_3 = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Verification/Ex_Vivo/2D/20180412_Unmasc/', num2str(i+2)];
    output_filename    = ['/home/', str_user, '/Dicom/DIXON/Verification/Ex_Vivo/2D/20180412_Unmasc/',num2str(i),'_',num2str(i+2),'.h5'];


    %% reading bruker acqp, method and fid files.

    ex1 = read_bru_experiment(acquisition_path_1);
    ex2 = read_bru_experiment(acquisition_path_2);
    ex3 = read_bru_experiment(acquisition_path_3);

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


    [ nX, nY, nZ ] = get_dimensions( ex1 );

    [ nX_acq, nY_acq, nZ_acq ] = get_dimensions_acq( ex1 );

    [ readout, E1, E2 ] = get_encoding_size( ex1 , nZ );

    [header] = fill_the_flexible_xml_header_special_3Flash(ex1, ex2, ex3);

    %% reshape the fid to match the fixed data structure and remove the zero

    [ data_for_acqp_1 ] = remove_zero_from_fid( ex1 );
    [ data_for_acqp_2 ] = remove_zero_from_fid( ex2 );
    [ data_for_acqp_3 ] = remove_zero_from_fid( ex3 );

    check_size(size(data_for_acqp_1),size(data_for_acqp_2),size(data_for_acqp_3));

    number_of_channels=size(data_for_acqp_1,2);
    number_of_lines=size(data_for_acqp_1,3);

    %% fill the fixed data structure

    [ idx1, flag1 ] = fill_the_idx( header , ex1);
    [ idx2, flag2 ] = fill_the_idx( header , ex2);
    [ idx3, flag3 ] = fill_the_idx( header , ex3);

    clear  check_indice

    for i = 1:size(idx1.kspace_encode_step_1,2)
        if (idx1.contrast(i)==0)
            check_indice(idx1.kspace_encode_step_1(i)+1,idx1.kspace_encode_step_2(i)+1,1)=1;
        end
    end

    figure()
    subplot(2,1,1); imagesc(check_indice(:,:,1));

    %% parallel imaging option

    [ Ysamp_regular , Ysamp_ACS , acceleration_factor_y] = check_options_for_parallel_imaging( ex1 );



    total_number_of_lines=number_of_lines*3;

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
    acqblock = ismrmrd.Acquisition(total_number_of_lines);

    % Set the header elements that don't change
    acqblock.head.version(:) = 1;
    acqblock.head.number_of_samples(:) = readout;
    acqblock.head.center_sample(:) = floor(readout/2);
    acqblock.head.active_channels(:) = number_of_channels;
    acqblock.head.read_dir  = repmat([1 0 0]',[1 total_number_of_lines]);
    acqblock.head.phase_dir = repmat([0 1 0]',[1 total_number_of_lines]);
    acqblock.head.slice_dir = repmat([0 0 1]',[1 total_number_of_lines]);

    % Loop over the acquisitions, set the header, set the data and append


    %% echo 1

    for acqno = 1:number_of_lines

        % Set the header elements that change from acquisition to the next
        % c-style counting
        acqblock.head.scan_counter(acqno) =  acqno-1;
        acqblock.head.idx.kspace_encode_step_1(acqno) = idx1.kspace_encode_step_1(acqno);
        acqblock.head.idx.kspace_encode_step_2(acqno) = idx1.kspace_encode_step_2(acqno);
        acqblock.head.idx.repetition(acqno) = idx1.repetition(acqno);
        acqblock.head.idx.slice(acqno) = idx1.slice(acqno);
        acqblock.head.idx.set(acqno) = idx1.set(acqno);
        acqblock.head.idx.contrast(acqno) = 0; %idx1.contrast(acqno);
        acqblock.head.sample_time_us(acqno)=1/ex1.acqp.SW_h*1e6;

        % Set the flags
        acqblock.head.flagClearAll(acqno);

        %     str_msg=sprintf('count %d Ysamp %d        ', acqno,   idx.kspace_encode_step_1(acqno)    ); disp( str_msg);

        if (acceleration_factor_y>1)

            if ismember(idx1.kspace_encode_step_1(acqno),Ysamp_ACS-1)
                %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne

                if ismember(idx1.kspace_encode_step_1(acqno),Ysamp_regular)
                    % both calibration and part of the undersampled pattern
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
                    %                 disp('both calibration and part of the undersampled pattern');
                else
                    % in ACS block but not part of the regular undersampling
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
                    %                 disp('in ACS block but not part of the regular undersampling');
                end
            end
        end

        if (flag1.first_in_encoding_step1(acqno)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
        end

        if (flag1.last_in_encoding_step1(acqno)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
        end

        if (flag1.first_in_encoding_step2(acqno)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP2', acqno);
        end

        if (flag1.last_in_encoding_step2(acqno)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP2', acqno);
        end

        if (flag1.first_in_repetition(acqno)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
        end

        if (flag1.last_in_repetition(acqno)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
        end

        if (flag1.first_in_slice(acqno)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
        end

        if (flag1.last_in_slice(acqno)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
        end

        % fill the data
        acqblock.data{acqno} = squeeze(data_for_acqp_1(:,:,acqno));

    end




    %% echo 2 

    for acqno = number_of_lines+1:number_of_lines*2

        index=acqno-number_of_lines;

        % Set the header elements that change from acquisition to the next
        % c-style counting
        acqblock.head.scan_counter(acqno) =  acqno-1;
        acqblock.head.idx.kspace_encode_step_1(acqno) = idx2.kspace_encode_step_1(index);
        acqblock.head.idx.kspace_encode_step_2(acqno) = idx2.kspace_encode_step_2(index);
        acqblock.head.idx.repetition(acqno) = idx2.repetition(index);
        acqblock.head.idx.slice(acqno) = idx2.slice(index);
        acqblock.head.idx.set(acqno) = idx2.set(index);
        acqblock.head.idx.contrast(acqno) = 1; % idx2.contrast(index);
        acqblock.head.sample_time_us(acqno)=1/ex2.acqp.SW_h*1e6;

        % Set the flags
        acqblock.head.flagClearAll(acqno);

        %     str_msg=sprintf('count %d Ysamp %d        ', acqno,   idx.kspace_encode_step_1(acqno)    ); disp( str_msg);

        if (acceleration_factor_y>1)

            if ismember(idx2.kspace_encode_step_1(index),Ysamp_ACS-1)
                %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne

                if ismember(idx2.kspace_encode_step_1(index),Ysamp_regular)
                    % both calibration and part of the undersampled pattern
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
                    %                 disp('both calibration and part of the undersampled pattern');
                else
                    % in ACS block but not part of the regular undersampling
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
                    %                 disp('in ACS block but not part of the regular undersampling');
                end
            end
        end

        if (flag2.first_in_encoding_step1(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
        end

        if (flag2.last_in_encoding_step1(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
        end

        if (flag2.first_in_encoding_step2(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP2', acqno);
        end

        if (flag2.last_in_encoding_step2(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP2', acqno);
        end

        if (flag2.first_in_repetition(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
        end

        if (flag2.last_in_repetition(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
        end

        if (flag2.first_in_slice(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
        end

        if (flag2.last_in_slice(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
        end

        % fill the data
        acqblock.data{acqno} = squeeze(data_for_acqp_2(:,:,index));

    end





    %% echo 3 

    for acqno = number_of_lines*2+1:number_of_lines*3

        index=acqno-number_of_lines*2;

        % Set the header elements that change from acquisition to the next
        % c-style counting
        acqblock.head.scan_counter(acqno) =  acqno-1;
        acqblock.head.idx.kspace_encode_step_1(acqno) = idx3.kspace_encode_step_1(index);
        acqblock.head.idx.kspace_encode_step_2(acqno) = idx3.kspace_encode_step_2(index);
        acqblock.head.idx.repetition(acqno) = idx3.repetition(index);
        acqblock.head.idx.slice(acqno) = idx3.slice(index);
        acqblock.head.idx.set(acqno) = idx3.set(index);
        acqblock.head.idx.contrast(acqno) = 2; % idx3.contrast(index);
        acqblock.head.sample_time_us(acqno)=1/ex3.acqp.SW_h*1e6;

        % Set the flags
        acqblock.head.flagClearAll(acqno);

        %     str_msg=sprintf('count %d Ysamp %d        ', acqno,   idx.kspace_encode_step_1(acqno)    ); disp( str_msg);

        if (acceleration_factor_y>1)

            if ismember(idx3.kspace_encode_step_1(index),Ysamp_ACS-1)
                %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne

                if ismember(idx3.kspace_encode_step_1(index),Ysamp_regular)
                    % both calibration and part of the undersampled pattern
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
                    %                 disp('both calibration and part of the undersampled pattern');
                else
                    % in ACS block but not part of the regular undersampling
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
                    %                 disp('in ACS block but not part of the regular undersampling');
                end
            end
        end

        if (flag2.first_in_encoding_step1(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
        end

        if (flag2.last_in_encoding_step1(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
        end

        if (flag2.first_in_encoding_step2(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP2', acqno);
        end

        if (flag2.last_in_encoding_step2(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP2', acqno);
        end

        if (flag2.first_in_repetition(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
        end

        if (flag2.last_in_repetition(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
        end

        if (flag2.first_in_slice(index)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
        end

        if (flag2.last_in_slice(index)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
        end

        % fill the data
        acqblock.data{acqno} = squeeze(data_for_acqp_3(:,:,index));

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

    % disp('sample time');
    % disp(1/ex.acqp.SW_h*1e6)

    close(figure(7))
    figure(7)
    for s = 1:4
        subplot(2,2,s)
        temp=squeeze(data_for_acqp_1(:,6,s:4:end));
        imagesc(abs(ifft_2D(temp))); colormap(gray);
    end
end



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