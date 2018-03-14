%% Read Bruker File

%clear all
%close all



[status,id]= system('whoami');

str_user= id(1:end-1);

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

%% Boucle reco

n_reco = [39];

for numb = 1:size(n_reco,2)
    % acquisition_path=['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Validation/RawData/In_Vitro/2D/No_Grappa/20171221/Dixon/24'];
    % output_tmp = ['/home/', str_user, '/Dicom/DIXON/Validation/RecoData/In_Vitro/2D/No_Grappa/20171221/Dixon/24'];
    %  acquisition_path=['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Validation/RawData/Ex_Vivo/2D/No_Grappa/20171221/Dixon/13'];
    %  output_tmp = ['/home/', str_user, '/Dicom/DIXON/Validation/RecoData/Ex_Vivo/2D/No_Grappa/20171221/Dixon/13'];
%     acquisition_path    = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Validation/RawData/Ex_Vivo/2D/No_Grappa/20180122/Coeur_CC/',num2str(n_reco(numb))];
%     output_tmp          = ['/home/', str_user, '/Dicom/DIXON/Validation/RecoData/Ex_Vivo/2D/No_Grappa/20180122/Coeur_CC/',num2str(n_reco(numb))];
%     acquisition_path    = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Validation/RawData/In_Vitro/2D/No_Grappa/20180227/',num2str(n_reco(numb))];
%     output_tmp          = ['/home/', str_user, '/Dicom/DIXON/Validation/RecoData/In_Vitro/2D/No_Grappa/20180227/',num2str(n_reco(numb))];
      acquisition_path    = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Validation/RawData/Ex_Vivo/2D/No_Grappa/20180216/',num2str(n_reco(numb))];
      output_tmp          = ['/home/', str_user, '/Dicom/DIXON/Validation/RecoData/Ex_Vivo/2D/No_Grappa/20180216/',num2str(n_reco(numb))];
%     acquisition_path    = ['/home/', str_user, '/mount/Imagerie/For_Kylian/Dixon/Validation/RawData/Dixon_t2star/Vitro/',num2str(n_reco(numb))];
%     output_tmp          = ['/home/', str_user, '/Dicom/DIXON/',num2str(n_reco(numb))];

    output_filename     = [output_tmp,'.h5'];


    %% reading bruker acqp, method and fid files. 

    ex = read_bru_experiment(acquisition_path);

    [ nX, nY, nZ ] = get_dimensions( ex );

    [ nX_acq, nY_acq, nZ_acq ] = get_dimensions_acq( ex );

    [ readout, E1, E2 ] = get_encoding_size( ex , nZ );

    [header] = fill_the_flexible_xml_header(ex);

    %% reshape the fid to match the fixed data structure and remove the zero

    [ data_for_acqp ]   = remove_zero_from_fid( ex );
    number_of_channels  = size(data_for_acqp,2);
    number_of_lines     = size(data_for_acqp,3);

    %% gestion de la diffusion

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

    [ idx, flag ] = fill_the_idx( header , ex);


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

    TE = ex.acqp.ACQ_echo_time;

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

        %     str_msg=sprintf('count %d Ysamp %d        ', acqno,   idx.kspace_encode_step_1(acqno)    ); disp( str_msg);

        if (acceleration_factor_y>1)

            if ismember(idx.kspace_encode_step_1(acqno),Ysamp_ACS-1)
                %% attention si on met Ysamp_ACS sans le moins 1 la phase n'est pas bonne

                if ismember(idx.kspace_encode_step_1(acqno),Ysamp_regular)
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

        if (flag.first_in_encoding_step1(acqno)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);

        end

        if (flag.last_in_encoding_step1(acqno)== 1)
            acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);

        end

         if (flag.first_in_contrast(acqno)== 1)
            acqblock.head.flagSet('ACQ_FIRST_IN_CONTRAST', acqno);
         end

          if (flag.last_in_contrast(acqno)== 1)
             acqblock.head.flagSet('ACQ_LAST_IN_CONTRAST', acqno);
          end
    %     
    %     

        %         acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
        %         acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
        %     elseif acqno==size(K,2)
        %         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
        %         acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
        %         acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
        %     end
    end
   
%      % Apply the offsets
         ch = size(data_for_acqp, 2);

        Offcenter = phaseOffsetCorrection(ex);

            data_for_acqp = permute(data_for_acqp,[3 1 2]);
            tmp = data_for_acqp;
        for ne=1:ex.method.PVM_NEchoImages
             for k=1:ch
                for jj=1:readout
                        data_for_acqp(ne:ex.method.PVM_NEchoImages:end,jj,k)= tmp(ne:ex.method.PVM_NEchoImages:end,jj,k).*Offcenter(:,jj,:);
                end
             end
        end
            data_for_acqp = permute(data_for_acqp,[2 3 1]);
         
            
   for acqno = 1:number_of_lines         
         acqblock.data{acqno} = squeeze(data_for_acqp(:,:,acqno));
   end
        
% Append the acquisition block
    dset.appendAcquisition(acqblock);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Fill the xml header %
    %%%%%%%%%%%%%%%%%%%%%%%%

    header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = acceleration_factor_y ;
    header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
    header.encoding.parallelImaging.calibrationMode = 'embedded' ;
    header.sequenceParameters.TE = TE;
    
    %% Serialize and write to the data set
    xmlstring = ismrmrd.xml.serialize(header);
    dset.writexml(xmlstring);

    %% Write the dataset
    dset.close();

    fprintf('-------------------------------------------------------------- \n');
    disp('fin de la conversion');
    fprintf('-------------------------------------------------------------- \n');
end