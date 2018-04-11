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



[status,id]= system('whoami');

str_user= id(1:end-1);

check_if_iam_using_the_ihuserver(str_user);

[ str_network_imagerie, str_network_perso ] = get_network_name( str_user );

% acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_1/4/';  % grappa 2 * 2
% acquisition_path='/home/valery/DICOM/badres/11/'; % fully  % ne marche
% pas ?
%  acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_1/15/'; % CS
% acquisition_path=['/home/',str_user, str_network_imagerie,'/DICOM_DATA/2017-07-07-test_mp2rage_1/16/'];  % fully
%  acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_1/19/'  % grappa 2 * 2
% acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_2/13'; %fully  ok
% acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_2/14'; %grappa 2  ok
%  acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_2/15'; %grappa 2   ok
% acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_2/16'; %grappa 2 suivant z  ok
% acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_2/17'; %grappa 2 * 2  
% acquisition_path='/home/valery/Reseau/Imagerie/DICOM_DATA/2017-07-07-test_mp2rage_2/18'; %CS 2 x 2
% acquisition_path='/home/valery/DICOM/Bruker_7617/6/';
% acquisition_path='/home/valery/DICOM/tests13062017/74/';
% acquisition_path='/home/valery/Reseau/Imagerie/Auckland/Kadence/Control/2016_05_03_Heart_2/17/';

% acquisition_path='/home/valery/DICOM/2017-10-23_mp2rage/2'; % phantom grappa
% acquisition_path='/home/valery/DICOM/2017-10-24_mp2rage_coeur/6';
% acquisition_path='/home/valery/DICOM/2017-10-24_mp2rage_coeur/7';
% acquisition_path='/home/valery/DICOM/99/'   %% human heart  grappa y

% acquisition_path='/home/valery/DICOM/103/'   %% human heart
acquisition_path='/home/valery/DICOM/105/'   %% human heart
% acquisition_path='/home/valery/DICOM/2017-10-24_mp2rage_coeur/8';
% acquisition_path='/home/valery/DICOM/2017-10-24_mp2rage_coeur/9';
output_filename = '/home/valery/DICOM/testdata.h5';

%% reading bruker acqp, method and fid files.

ex = read_bru_experiment(acquisition_path);

[ nX, nY, nZ ] = get_dimensions( ex );

[ readout, E1, E2 ] = get_encoding_size( ex , nZ )

[header] = fill_the_flexible_xml_header(ex);


%% reshape the fid to match the fixed data structure


[readout, number_of_channels]=get_additionnal_info_for_fid(ex);

check_if_zero_problem( ex );

ADC_size=128;

sizeR2=ADC_size*ceil(readout*number_of_channels/ADC_size);

acqp_begin=reshape(ex.fid, sizeR2,[]);

acqp_end=acqp_begin(1:(readout*number_of_channels),:);

number_of_lines=size(acqp_end,2);

data_for_acqp=reshape(acqp_end,[readout,number_of_channels,number_of_lines]);

number_of_echos=ex.method.PVM_NEchoImages;

if (size(data_for_acqp,3)/number_of_echos~=ex.method.Genpoint)
    error('problem dans le nombre de ligne')
end

if (strcmp(ex.method.EnableCS,'Fully') ||  strcmp(ex.method.EnableCS,'No'))
    
    disp('Fully')
    
    ex.method.AccCSY=1;
    
    ex.method.AccCSZ=1;
    
    acceleration_factor_y=ex.method.AccCSY;
    
    acceleration_factor_z=ex.method.AccCSZ;
    
elseif (strcmp(ex.method.EnableCS,'Grappa'))
    
    disp('Grappa')
    
    acceleration_factor_y=ex.method.AccCSY;
    
    acceleration_factor_z=ex.method.AccCSZ;
    
elseif (strcmp(ex.method.EnableCS,'Yes') || strcmp(ex.method.EnableCS,'Compressed_Sensing'))  % Compress_Sensing
    
    disp('Compress Sensing')
    acceleration_factor_y=ex.method.AccCSY;
    
    acceleration_factor_z=ex.method.AccCSZ;
else
    
    error('problem option');
    
end





%% fill the fixed data structure

[ idx ] = fill_the_idx_mp2rage( header , ex);

if (mod(readout,2)==1)
    
    readout=readout+1;
    header.encoding.encodedSpace.matrixSize.x=readout;
    header.encoding.reconSpace.matrixSize.x=readout;
    
    kspace=zeros(readout, size( data_for_acqp,2),size( data_for_acqp,3)); 
    
    kspace(1:end-1,:,:)=data_for_acqp(:,:,:);     
      
    clear data_for_acqp
   
    data_for_acqp=kspace;
    
    clear  kspace
    
end



clear  check_indice

for i = 1:size(idx.kspace_encode_step_1,2)
    if (idx.contrast(i)==0)
        check_indice(idx.kspace_encode_step_1(i)+1,idx.kspace_encode_step_2(i)+1,1)=1;
    end
    
    if (idx.contrast(i)==1)
        check_indice(idx.kspace_encode_step_1(i)+1,idx.kspace_encode_step_2(i)+1,2)=1;
    end
end


figure()
subplot(2,1,1); imagesc(check_indice(:,:,1));
subplot(2,1,2); imagesc(check_indice(:,:,2));


% [ idx ] = fill_the_idx_mp2rage_test( header , ex);

%% parallel imaging option

% [ Ysamp_regular , Ysamp_ACS ] = check_options_for_parallel_imaging( ex );


if (acceleration_factor_y>1 || acceleration_factor_z>1)
    
    
    ligne_du_milieu_suivant_y=check_indice(:,round(size(check_indice,2)/2),1);
    
    clear ligne_du_milieu_suivant_y_value
    
    for j = 1:size(ligne_du_milieu_suivant_y,1)
        
        if (ligne_du_milieu_suivant_y(j)==1)
            ligne_du_milieu_suivant_y_value(j)=j;
        end
        
    end
    
    ligne_du_milieu_suivant_y_value(ligne_du_milieu_suivant_y_value==0)=[];
    
    
    
    ligne_du_milieu_suivant_z=check_indice(round(size(check_indice,1)/2),:,1);
    
    clear ligne_du_milieu_suivant_z_value
    
    for j = 1:size(ligne_du_milieu_suivant_z,2)
        
        if (ligne_du_milieu_suivant_z(j)==1)
            ligne_du_milieu_suivant_z_value(j)=j;
        end
        
    end
    
    ligne_du_milieu_suivant_z_value(ligne_du_milieu_suivant_z_value==0)=[];
    
    % on recupere les numero de lignes acs suivants y
    if (acceleration_factor_y>1 && acceleration_factor_z==1)
        
        tempo=diff(ligne_du_milieu_suivant_y);
        acs_indice_y=find(tempo==0);
        regular_sampling_y=[1:2:nY];
    end
    
    if (acceleration_factor_y==1 && acceleration_factor_z>1)
        
        tempo=diff(ligne_du_milieu_suivant_z);
        acs_indice_z=find(tempo==0);
        regular_sampling_z=[1:2:nZ];
        
    end
    
     if (acceleration_factor_y>1 && acceleration_factor_z>1)
         
%           tempo=diff(ligne_du_milieu_suivant_y);
        acs_indice_y=find(ligne_du_milieu_suivant_y==1);
        regular_sampling_y=[1:2:nY];
         
%          tempo=diff(ligne_du_milieu_suivant_z);
        acs_indice_z=find(ligne_du_milieu_suivant_z==1);
        regular_sampling_z=[1:2:nZ]; 
        
        
     end
    
    
%     
%     
%     ACShw=round(ex.method.CenterMaskSize/2);
%     
%     clear  Ysamp_u Ysamp_ACS Ysamp
%     
%     Ysamp_u = [1:2:nY] ; % undersampling by every alternate line
%     Ysamp_regular=Ysamp_u;
%     Ysamp_ACS = [nY/2-ACShw+1 : nY/2+ACShw] ; % GRAPPA autocalibration lines
%     Ysamp = union(Ysamp_u, Ysamp_ACS) ; % actually sampled lines
%     nYsamp = length(Ysamp) ; % number of actually sampled
    
    %     figure(2)
    %     plot(Ysamp, 'r');  hold on ; plot(ligne_du_milieu_suivant_y_value+0.5, 'b')
    
    
    
end
% 
% kspace=zeros(200,7,200,200,2);
% 
% 
% for acqno = 1:number_of_lines
%     
%         
%    e1=idx.kspace_encode_step_1(acqno)+1;
%    e2=idx.kspace_encode_step_2(acqno)+1;
%    echo=idx.contrast(acqno)+1;    
%    
%    kspace(:,:,e1,e2, echo)=   data_for_acqp(:,:,acqno); 
%    
%    
%    
% end
% 
% 
% for acqno = 1:number_of_lines
%     
%     e1=idx.kspace_encode_step_1(acqno)+1;
%     e2=idx.kspace_encode_step_2(acqno)+1;
%     echo=idx.contrast(acqno)+1;    
%     
%     new_data_for_acqp(:,:,acqno)= kspace(:,:,e1,e2, echo);
% 
% end



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
    acqblock.head.idx.contrast(acqno) = idx.contrast(acqno);
    acqblock.head.sample_time_us(acqno)=1;    
    %% TODO a verifier EffReadOffset ?
    acqblock.head.position(:,acqno)=[ex.method.PVM_EffPhase0Offset,ex.method.PVM_EffPhase1Offset,ex.method.PVM_EffPhase2Offset];
    
    % Set the flags
    acqblock.head.flagClearAll(acqno);
    
    %     str_msg=sprintf('count %d Ysamp %d        ', acqno,   idx.kspace_encode_step_1(acqno)    ); disp( str_msg);
    
    if (acceleration_factor_y>1 && acceleration_factor_z==1)
        
        if ismember(idx.kspace_encode_step_2(acqno),ligne_du_milieu_suivant_z_value-1)
            
            if ismember(idx.kspace_encode_step_1(acqno),acs_indice_y-1)
                %                 if ismember(idx.kspace_encode_step_1(acqno),Ysamp_ACS-1)
                
                
                %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne
                %la conditon dans le carrée est respecté
                
                if ismember(idx.kspace_encode_step_1(acqno),regular_sampling_y-1)
                    % both calibration and part of the undersampled pattern
                    str_msg=sprintf('y> 1 & z==1  %d ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING %d %d', acqno , idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_msg);
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
                    str_smg=sprintf('y>1 && z==1  ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING %d  %d  %d', acqno, idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_smg);
                    %                 disp('both calibration and part of the undersampled pattern');
                else
                    % in ACS block but not part of the regular undersampling
                    str_msg=sprintf('y> 1 & z==1  %d ACQ_IS_PARALLEL_CALIBRATION %d %d', acqno , idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_msg);
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
                    str_smg=sprintf('y>1 && z==1  ACQ_IS_PARALLEL_CALIBRATION %d  %d  %d ', acqno, idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_smg);
                    %                 disp('in ACS block but not part of the regular undersampling');
                end
                 
            end
        end
    end
    
    
    if (acceleration_factor_y==1 && acceleration_factor_z>1)
        
        if ismember(idx.kspace_encode_step_1(acqno),ligne_du_milieu_suivant_y_value-1)
            
            if ismember(idx.kspace_encode_step_2(acqno),acs_indice_z-1)
                %                 if ismember(idx.kspace_encode_step_1(acqno),Ysamp_ACS-1)                
                
                %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne
                %la conditon dans le carrée est respecté
                
                if ismember(idx.kspace_encode_step_2(acqno),regular_sampling_z-1)
                    % both calibration and part of the undersampled pattern
                    str_msg=sprintf('y==1 & z>1  %d ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING %d %d', acqno , idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_msg);
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
                    str_smg=sprintf('y==1 && z>1  ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING %d  %d  %d ', acqno, idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_smg);
                    %                 disp('both calibration and part of the undersampled pattern');
                else
                    % in ACS block but not part of the regular undersampling
                    str_msg=sprintf('y==1 & z>1  %d ACQ_IS_PARALLEL_CALIBRATION %d %d', acqno , idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_msg);
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
                    str_smg=sprintf('y==1 && z>1  ACQ_IS_PARALLEL_CALIBRATION %d  %d  %d ', acqno, idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_smg);
                    %                 disp('in ACS block but not part of the regular undersampling');
                end
                
                
                
            end
        end
    end
    
    
    if (acceleration_factor_y>1 && acceleration_factor_z>1)
        
        if ismember(idx.kspace_encode_step_1(acqno),acs_indice_y-1)
            
            if ismember(idx.kspace_encode_step_2(acqno),acs_indice_z-1)
                %                 if ismember(idx.kspace_encode_step_1(acqno),Ysamp_ACS-1)                
                
                %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne
                %la conditon dans le carrée est respecté
                
                if ismember(idx.kspace_encode_step_2(acqno),regular_sampling_z-1) && ismember(idx.kspace_encode_step_1(acqno),regular_sampling_y-1)
                    % both calibration and part of the undersampled pattern
                     str_msg=sprintf('y> 1 & z>1  %d ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING %d %d', acqno , idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_msg);
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
                    %                 disp('both calibration and part of the undersampled pattern');
                else
                    % in ACS block but not part of the regular undersampling
                     str_msg=sprintf('y> 1 & z>1  %d ACQ_IS_PARALLEL_CALIBRATION %d %d', acqno , idx.kspace_encode_step_1(acqno), idx.kspace_encode_step_2(acqno)); disp(str_msg);
                    acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
                    %                 disp('in ACS block but not part of the regular undersampling');
                end  
                
            end
        end
    end
    
    
    
    %     if (acceleration_factor_z>1)
    %
    %         if ismember(idx.kspace_encode_step_2(acqno),ligne_du_milieu_suivant_z_value-1)
    %
    %             if ismember(idx.kspace_encode_step_1(acqno),Ysamp_ACS-1)
    %                 %% attention si on mais Ysamp_ACS sans le moins 1 la phase n'est pas bonne
    %
    %                 if ismember(idx.kspace_encode_step_1(acqno),Ysamp_regular-1)
    %                     % both calibration and part of the undersampled pattern
    %                     acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', acqno);
    %                     %                 disp('both calibration and part of the undersampled pattern');
    %                 else
    %                     % in ACS block but not part of the regular undersampling
    %                     acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', acqno) ;
    %                     %                 disp('in ACS block but not part of the regular undersampling');
    %                 end
    %             end
    %         end
    %     end
    
    
    
    %     if (flag.first_in_encoding_step1(acqno)== 1)
    %         acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', acqno);
    %     end
    %
    %     if (flag.last_in_encoding_step1(acqno)== 1)
    %         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
    %     end
    
    %         acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
    %         acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
    %     elseif acqno==size(K,2)
    %         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', acqno);
    %         acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
    %         acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
    %     end
    
    % fill the data
    acqblock.data{acqno} = squeeze(data_for_acqp(:,:,acqno));
    
end

% Append the acquisition block
dset.appendAcquisition(acqblock);



%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%

header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = acceleration_factor_y ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = acceleration_factor_z ;
header.encoding.parallelImaging.calibrationMode = 'embedded' ;


%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();



disp('sample time');
disp(1/ex.acqp.SW_h*1e6)





%%%%%%%%%%%%
%% DEBUG
% rare_factor=ex.acqp.ACQ_rare_factor;
%
% test_reco=reshape(data_for_acqp,readout,number_of_channels,rare_factor,number_of_echos,[]);
% test_reco=permute(test_reco,[1,2,3,5,4]);
% test_reco=reshape(test_reco,readout,number_of_channels,[],number_of_echos);
%
% test_reco=reshape(data_for_acqp,readout,number_of_channels,E1,number_of_echos,E2);
% % test_reco=permute(test_reco,[1,3,5,2,4]);
%
% size(test_reco)
%
% data_for_acqp_test=reshape( test_reco,  readout, number_of_channels, []);

% data_for_acqp_test2=zeros(size(data_for_acqp_test,1)+1,size(data_for_acqp_test,2),size(data_for_acqp_test,3));

% data_for_acqp_test2(1:end-1,:,:)=data_for_acqp_test;


% echo1=ifft_3D(test_reco(:,:,:,:,1));
% echo1_rss=sqrt(sum(abs(echo1).^2,4));
%
% imagesc(echo1_rss(:,:,50))


%%%%%%%%%%%%


% size(test_reco)

%
% GradPhaseVector=ex.method.GradPhaseVector*E1/2;
% GradPhaseVector=round(GradPhaseVector-min(GradPhaseVector)+1);
%
% GradSliceVector=ex.method.GradSliceVector*E2/2;
% GradSliceVector=round(GradSliceVector-min(GradSliceVector)+1);
%
% count=0;
%
% for  ne = 1:1:number_of_echos
%     for  ll = 1:1:size(GradPhaseVector,1)
%
%         count=count+1;
%
%         str_msg=sprintf('count %d e1 %d e2 %d   ', count, GradPhaseVector(ll), GradSliceVector(ll)); disp( str_msg);
%
%         idx.kspace_encode_step_1(count)=GradPhaseVector(ll);
%         idx.kspace_encode_step_2(count)=GradSliceVector(ll);
%
%         idx.slice(count)=0;
%         idx.repetition(count)=0;
%         idx.contrast(count)=ne-1;
%
%     end;
% end

%rare factor

% plot(abs(ex.fid(1:10000)))
%
% count=0;
%
% matrice=zeros(readout,E1,E2,number_of_channels,number_of_echos);
%
% for  e2 = 1:1:E2
%     for  ne = 1:1:number_of_echos
%         for  e1 = 1:1:E1
%             count=count+1;
% %             for  c = 1:1:number_of_channels
%
%
%             matrice(:,e1,e2,:,ne)=data_for_acqp(:,:,count);
% %             end
%         end
%     end
% end
%
%
% imagesc(abs(matrice(:,:,50,1,1)));
%
% echo1=ifft_3D(matrice(:,:,:,:,1));
%
% echo2_rss=sqrt(sum(abs(echo1).^2,4));
