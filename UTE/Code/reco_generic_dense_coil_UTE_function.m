function [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  lili)


%% NOTATION Matlab Gadgetron à completer


% - number_of_output_respiratory_phases
% - number_of_output_cardiac_phases
% - trajectory_subsets
% - kspace_subsets
% - profiles_per_frame
% - calculate_trajectory_for_reconstruction
% - calculate_density_compensation_for_reconstruction
% - calculate_trajectory_for_frame
% - number_of_projections
% - get_maximum_umber_of_projections_index
% - set_num_projections_per_batch
% - projections_percentage
% - num_projections_expected_
% - num_projections_to_use_
% - projection_subsets
% - projections_per_recon_

%% TODO

% la detection de peak cardiaque peut être améliorée
% ajouter la MOCO


% clear all
% close all

addpath('../Generic_functions/');

[ str_user ] = get_PC_name();

% folder= ['/home/', str_user ,
% '/Reseau/Imagerie/DICOM_DATA/2017-03-08-UTE/'];
folder= ['/home/', str_user , '/DICOM/2017-03-08-UTE/'];
name = 'meas_MID00393_FID16077_UTE_SG1_5_postinf_FA30_NR10_2.h5';

filename = [ folder, '/FID/'  name ] ;

if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end


folder_data='/tmp/gadgetron/';


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
    nbSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    nbSlices = 1;
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
    nPhases = hdr.encoding.encodingLimits.phase.maximum + 1;
catch
    nPhases = 1;
end

try
    nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 + 1;
catch
    nContrasts = 1;
end




%% quelques vérifications

if (hdr.encoding.encodedSpace.matrixSize.x~=hdr.encoding.reconSpace.matrixSize.x)
    
    disp(' presence d oversampling \n');
end


%% partie reco
ReadOffCenter  = hdr.userParameters.userParameterDouble(1,5).value; % input for phase roll
PhaseOffCenter = hdr.userParameters.userParameterDouble(1,6).value; % input for phase roll
acqOffSet      = hdr.userParameters.userParameterDouble(1,7).value; % input for number of k=0 point (selfgating)
AcqMode        = hdr.userParameters.userParameterDouble(1,8).value; % input : NONE ,REORDER, GODLEN ANGLE

%% partie physio
info.study_time=hdr.studyInformation.studyTime;
info.study_time_read= strrep(info.study_time,':','');
info.Tr_=hdr.sequenceParameters.TR;
nCh=hdr.acquisitionSystemInformation.receiverChannels;
nbProjections=hdr.encoding.encodingLimits.kspace_encoding_step_1.maximum -hdr.encoding.encodingLimits.kspace_encoding_step_1.minimum +1 ;


%% parametres

% frequence du signal ute
Fs_ute = 1/(info.Tr_/1000); % in Hz

number_of_projections=nbProjections;
readout=enc_Nx;

number_of_channels=30;
channels=8;
repetition=10;


%% parametres du self gating , parametre qui doit venir de la special card 

threshold.exclusion=1000;
threshold.debut=threshold.exclusion+1;
threshold.fin=number_of_projections-threshold.exclusion;

threshold.resp=1;
threshold.ecg=1;
threshold.svd=0.1;

threshold.exclusion_all=1000;
threshold.debut_all=threshold.exclusion_all+1;
threshold.fin_all=number_of_projections*repetition-threshold.exclusion_all;

threshold.exclusion_ecg=1000;
threshold.debut_ecg=threshold.exclusion_ecg+1;
threshold.fin_ecg=number_of_projections*repetition-threshold.exclusion_ecg/2;


%% chargement des données de self gating

data.input.raw=zeros(number_of_projections,number_of_channels,repetition);

for r=1:1:repetition
    
    str_r=sprintf('%d',r-1);
    
    filename=[folder_data,'/dense_data_store_abs_',str_r,'.dat'];
    disp(filename);
    temp=load(filename);
    data.input.raw(:,:,r)=temp;
    
end


%% on tri les données et on applique le filtre sur toutes les projections

data.allrep.raw=reformat_to_one_block(data.input.raw);

data.allrep.filter.resp=apply_butterworth_filtering(data.allrep.raw,1, Fs_ute);

data.allrep.substract.resp=zeros(size(data.allrep.raw));

for c=1:1:number_of_channels
    data.allrep.substract.resp(:,c)=data.allrep.raw(:,c)-data.allrep.filter.resp(:,c);
end

%% on reformate les données par répétition

[ data.filter.resp ] = reformat_to_each_repitition( data.allrep.filter.resp, repetition );

[ data.substract.resp ] = reformat_to_each_repitition( data.allrep.substract.resp, repetition );

%% ici on boucle sur les repetitions

[ navigator ] = extract_navigator_resp_and_ecg( data , threshold.debut, threshold.fin ,repetition , threshold , Fs_ute );

navigator.allrep.resp.normalize=single(navigator.resp.normalize(:));

navigator.allrep.ecg.filter=single(navigator.ecg.filter(:));


%% partie respiration

 % parametre qui doit venir de la special card 

[ table ] = extract_liste_respiratory( single(navigator.allrep.resp.normalize), threshold.debut_all, threshold.fin_all, number_of_output_respiratory_phases );


%% partie ecg

%% on prend manuellement les coils 21 et 22 car sur les repetitions 8 et 9 l'auto dectection ne marche pas

navigator.ecg.test=  (data.allrep.substract.resp(:,21)+ data.allrep.substract.resp(:,22))/2;

[ navigator.ecg.testfilter(:,1) ]=apply_butterworth_filtering(navigator.ecg.test(:,1),3, Fs_ute);



%% detection des peaks

[ data_peak_dectection.plus, data_peak_dectection.moins ] = find_peaks(single(navigator.ecg.testfilter), threshold.debut_ecg , threshold.fin_ecg);

[ time_difference.plus ] = find_time_difference_between_peaks( data_peak_dectection.plus, info.Tr_ );

[ time_difference.moins ] = find_time_difference_between_peaks( data_peak_dectection.moins, info.Tr_ );

% ca c'est just pour enlever le delai entre les pics + et -
delay=data_peak_dectection.plus(:,1)-data_peak_dectection.moins(:,1);

data_peak_dectection.plusdelay(:,1)=data_peak_dectection.plus(:,1)-median(delay);
data_peak_dectection.plusdelay(:,2)=data_peak_dectection.plus(:,2);

%% on enleve les peaks qui ont mal été détecté


threshold.repolarisation_time=300;
threshold.maximum_time=800;


data_peak_dectection.plus_new=remove_bad_peaks(data_peak_dectection.plus , time_difference.plus , threshold.repolarisation_time, info.Tr_ , threshold.maximum_time);

[ time_difference.plus_new ] = find_time_difference_between_peaks( data_peak_dectection.plus_new, info.Tr_ );

data_peak_dectection.moins_new=remove_bad_peaks(data_peak_dectection.moins , time_difference.moins , threshold.repolarisation_time, info.Tr_ , threshold.maximum_time);

[ time_difference.moins_new ] = find_time_difference_between_peaks( data_peak_dectection.moins_new, info.Tr_ );

%% nombre de répétions qui seront utilisées
% on peut diminuer ce chiffre si on manque de mémoire

repetition_for_reco=10;

%% sort ecg data

%g.NumberofPhase  % parametre qui doit venir de la special card 

[ table_ecg_self_gating ] = extract_list_ecg( single(data_peak_dectection.moins_new(:,1)), number_of_output_cardiac_phases ,repetition_for_reco, number_of_projections  );


%% nouvelle numerotation en partant de la systole

debut_change=5;

[ table_ecg_self_gating_in_systole ] = start_list_ecg_in_systole( table_ecg_self_gating , debut_change, number_of_output_cardiac_phases);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% lecture de ECG enregistré
% parametre qui doit venir de la special card 
% filename=['/home/', str_user , '/Reseau/Imagerie/DICOM_DATA/2017-03-08-UTE/ecg/UTEPMUsignal_20170309T',info.study_time_read,'.mat'];
filename=['/home/', str_user , '/DICOM/2017-03-08-UTE/ecg/UTEPMUsignal_20170309T',info.study_time_read,'.mat'];

[ ecg_truth ] = read_ecg_from_mat( filename );

T_ecg = 1/400;            % Sample time
L_ecg = size(ecg_truth,1);  % Length of signal
t_ecg = (0:L_ecg-1)*T_ecg;        % Time vector


%% calcul de vrai Tr de l'ute d'après l'ECG
%t_ecg(end)/300000*1000

T_ute = t_ecg(end)/(number_of_projections*10)   ; %1167/300000*1000; % Sample time
L_ute = size(navigator.ecg.testfilter,1);  % Length of signal
t_ute = (0:L_ute-1)*T_ute;        % Time vector

%% resampling of data peak detection

peak_detection.plus_true_timing(:,1)=data_peak_dectection.plus(:,1)*T_ute;
detection.moins_true_timing(:,1)=data_peak_dectection.moins(:,1)*T_ute;

peak_detection.plus_true_timing(:,2)=data_peak_dectection.plus(:,2);
peak_detection.moins_true_timing(:,2)=data_peak_dectection.moins(:,2);

%% resampling of the recorded ecg in the ute freqency acquisition

[ ecg_resample ] = resample_ecg_to_ute_frequency( ecg_truth , number_of_projections*10);

%% peak detection

[pks,locs]=findpeaks(-double(ecg_resample(:,2)),'MinPeakDistance',10,'MinPeakHeight', 200 );
clear pks

[ table_ecg_recording ] = extract_list_ecg( locs(:,1), number_of_output_cardiac_phases ,repetition_for_reco, number_of_projections  );


%% temporal filter
%% parametre qui doivent venir de la special card 

mode.filter=lili;

mode.filter

if (strcmp(mode.filter,'Y'))
    temporal_filter.index=3;
else
    temporal_filter.index=1;
end

temporal_filter.vecteur_ratio=[3 2 1]/4;

temporal_filter.size=size(temporal_filter.vecteur_ratio,2);

[ table_ecg_recording_extended ] = extract_list_ecg( locs(:,1), number_of_output_cardiac_phases*temporal_filter.index ,repetition_for_reco, number_of_projections  );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save physio for gnuplot

% save_physio_for_gnuplot( t_ecg(1,:)', ecg_truth(:,2) * -1/50  , ['/home/', str_user , '/Reseau/Valery/MatlabUnix/GPU/MakeFigures/'], 'ecg_voie1.txt' );
save_physio_for_gnuplot( t_ecg(1,:)', ecg_truth(:,2) * -1/50  , ['/tmp/gadgetron/'], 'ecg_voie1.txt' );
% save_physio_for_gnuplot( t_ecg(1,:)', ecg_truth(:,3) * 1/50  , ['/home/', str_user , '/Reseau/Valery/MatlabUnix/GPU/MakeFigures/'], 'ecg_voie2.txt' );
save_physio_for_gnuplot( t_ecg(1,:)', ecg_truth(:,3) * 1/50  , ['/tmp/gadgetron/'], 'ecg_voie2.txt' );

%% plot the respiratory navigator
plot_figure_resp( table,navigator, 1, 10000 );

%% plot the ecg navigator
plot_figure_ecg( t_ecg , ecg_truth, t_ute , navigator, data_peak_dectection,T_ute, 20, 40 );

%% plot the ecg navigator with the cardiac phase
plot_figure_ecg_sort_phase( t_ecg, ecg_truth, t_ute , navigator , data_peak_dectection, T_ute, table_ecg_recording, table_ecg_self_gating_in_systole ,  20, 40);

%% plot the ecg recording
plot_figure_ecg_recording( t_ute, ecg_resample, locs , T_ute, 20, 40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% lecture du kspace apres coil compression

% TODO chnager le nom de data input
% TODO afficher la liste avec le nombre de projection par repetition
repetition_for_reco=10;

input.matrice=zeros(readout, number_of_projections ,channels, repetition_for_reco);

for r=1:1:repetition_for_reco
    
    t = cputime;
    
    str_r=sprintf('%d',r-1);
    
    filename= ['/home/', str_user ,'/Tempo/savekspace_ute_',str_r,'.bin']
    [ tempo ] =single( read_binary_complex( filename, [ readout*channels* number_of_projections ] ));
    
    tempo1=reshape(tempo, [readout, channels,  number_of_projections]);
    
    input.matrice(:,:,:,r)= permute(tempo1, [1 3 2]);
    
    e = cputime-t;
    str_msg=sprintf('%d %f  %f %f \n', r, e , (e -mod(e,60))/60 , mod(e,60)   ); disp(str_msg);
    
end

clear tempo tempo1 data

disp('fin lecture des données');


%% parametre qui doivent venir de la special card

mode.recording='Y';
mode.selfgating='N';

if (strcmp(mode.selfgating,'Y'))
    
    table_ecg_for_reco=table_ecg_self_gating_in_systole;
    
elseif (strcmp(mode.recording,'Y'))
    
    table_ecg_for_reco=table_ecg_recording;
    
else
    
end


%% parametre qui doivent venir de la special card ou du header m1

info.N=2*size(input.matrice,1);
info.nSl=info.N;

info.sample_time_us = 5; %% TODO a verifier et à rendre générique
info.osf       = 1;  %oversampling
info.SGPoints  = 2*acqOffSet;
info.GA        = AcqMode;
info.RiseT     = 300;
info.BWp       = (1/(info.sample_time_us*1e-6))/(info.osf*info.N); % Hz/px



clear kspace_selecting trajectory_selecting 
clear trajectory_selecting_with_filter  trajectory_selecting_with_filter
clear FT mattest DCF matrice_nRep traj_nRep  temp ktrajIdxAvant ktrajIdxApres ktrajIdx image_nocrop image_crop



for repetition_for_reco_indice=10:1:repetition_for_reco;
    
    str_msg=sprintf(' reconstruction  %d / %d   \n', repetition_for_reco_indice ,  repetition_for_reco ); disp(str_msg);
    
    [enc_Nx,number_of_projections,nCh,nRep_use]=size(input.matrice(:,:,:,1:repetition_for_reco_indice));
    
    str_msg=sprintf(' nPE  %d nCh %d  nFE %d  nRep_use %d \n', enc_Nx,nCh,number_of_projections, nRep_use ); disp(str_msg);
    
    clear traj
    
    traj = ComputeUTEtraj(info.N, number_of_projections , info.SGPoints, info.osf,info.BWp,info.RiseT, info.GA);
    
    clear matrice_nRep  traj_nRep
    
    matrice_nRep=input.matrice(:,:,:,1:repetition_for_reco_indice);
    
    %% Gather repetition
    % data
    matrice_nRep=permute(matrice_nRep,[1 3 2 4]);
    
    matrice_nRep=matrice_nRep(:,:,:);
    
    matrice_nRep=permute(matrice_nRep,[1 3 2]);
    
    str_msg=sprintf(' nPE  %d nCh %d  nFE %d  nRep_use %d \n', enc_Nx,nCh,number_of_projections,nRep_use ); disp(str_msg);
    
    %% partie trajectoire
    traj_nRep=traj;
    
    for i=1:nRep_use-1
        traj_nRep=[traj_nRep; traj];
    end
    
    traj_nRep=reshape(traj_nRep,enc_Nx,[],3);
    
    clear traj
    
    str_msg=sprintf(' traj_nRep    %d   %d   %d \n',  size(traj_nRep,1),size(traj_nRep,2),size(traj_nRep,3) ); disp(str_msg);
    str_msg=sprintf(' matrice_nRep    %d   %d   %d \n',  size(matrice_nRep,1),size(matrice_nRep,2),size(matrice_nRep,3) ); disp(str_msg);
    
    %% SORT DATA
    
    kspace_subsets=cell(number_of_output_cardiac_phases*3,number_of_output_respiratory_phases);
    trajectory_subsets=cell(number_of_output_cardiac_phases*3,number_of_output_respiratory_phases);
    
    %% Image Reconstruction    
    
    for j=1:1:number_of_output_respiratory_phases
        
        for i=1:(number_of_output_cardiac_phases*temporal_filter.index)
            
            clear ktrajIdx
            %% on met a zero les projections qui ne correspondent pas à la bonne phase respiratoire
            [ table_input_ecg ] = remove_kspace_from_other_respiratory_stage(table.respiration.liste_ready(1:(number_of_projections*repetition_for_reco_indice)) , j,  table_ecg_recording_extended );
            ktrajIdx=find(table_input_ecg(:)==i);
            trajectory_subsets{i,j}=traj_nRep(:,ktrajIdx,:);
            kspace_subsets{i,j}=matrice_nRep(:,ktrajIdx,:);
            
        end
    end
    
    %% desallocation
    
    clear traj_nRep matrice_nRep
        
    %% if temporal filter yes
    if (strcmp(mode.filter,'Y'))
        
        
        %%  Realloc filt
        %je pars du principe qu'on a augment? le nombre de phaseECG par x3
        kspace_subsets_with_filter=cell(length(2:3:(number_of_output_cardiac_phases*temporal_filter.index)-1),number_of_output_respiratory_phases);
        trajectory_subsets_with_filter=cell(length(2:3:(number_of_output_cardiac_phases*temporal_filter.index)-1),number_of_output_respiratory_phases);
        
        for j=1:number_of_output_respiratory_phases
            count=0;
            for i=2:temporal_filter.index:(number_of_output_cardiac_phases*temporal_filter.index)-1
                % Ajout des 3 phases les unes ? c?t? des autres
                
                count=count+1;
                
                kspace_subsets_with_filter{count,j}=[reshape(kspace_subsets{i-1,j},[],8);reshape(kspace_subsets{i,j},[],8); reshape(kspace_subsets{i+1,j},[],8)];
                trajectory_subsets_with_filter{count,j}=[reshape(trajectory_subsets{i-1,j},[],3);reshape(trajectory_subsets{i,j},[],3); reshape(trajectory_subsets{i+1,j},[],3)];
                
                %Si on veut r?cup?rer une partie des donn?es des voisins.
                for ll=0:temporal_filter.size-1
                    %                 kspace_subsets_with_filter{count,j}=[kspace_subsets_with_filter{count,j};reshape(kspace_subsets{mod((i-1)-2,60)+1,j}(round(3/4*end):end,:),[],1);reshape(kspace_subsets{mod((i-1)+2,60)+1,j}(round(3/4*end):end,:),[],1)];
                    %                 trajectory_subsets_with_filter{count,j}=[trajectory_subsets_with_filter{count,j};reshape(trajectory_subsets{mod((i-1)-2,60)+1,j}(round(3/4*end):end,:,:),[],3);reshape(kspace_subsets{mod((i-1)+2,60)+1,j}(round(3/4*end):end,:,:),[],3)];
                    fprintf('%d Complet : %d + %d + %d et %f des hautes frequences de %d et %d\n',count, i-1,i,i+1,temporal_filter.vecteur_ratio(ll+1) ,...
                        mod((i-1)-ceil(temporal_filter.index/2)-ll,60)+1,mod((i-1)+ceil(temporal_filter.index/2)+ll,60)+1);
                    
                    kspace_subsets_with_filter{count,j}=[kspace_subsets_with_filter{count,j}; ...
                        reshape(kspace_subsets{mod((i-1)-ceil(temporal_filter.index/2)-ll,60)+1,j}(round(temporal_filter.vecteur_ratio(ll+1)*end):end,:),[],8);...
                        reshape(kspace_subsets{mod((i-1)+ceil(temporal_filter.index/2)+ll,60)+1,j}(round(temporal_filter.vecteur_ratio(ll+1)*end):end,:),[],8)];
                    
                    trajectory_subsets_with_filter{count,j}=[trajectory_subsets_with_filter{count,j};...
                        reshape(trajectory_subsets{mod((i-1)-ceil(temporal_filter.index/2)-ll,60)+1,j}(round(temporal_filter.vecteur_ratio(ll+1)*end):end,:),[],3);...
                        reshape(trajectory_subsets{mod((i-1)+ceil(temporal_filter.index/2)+ll,60)+1,j}(round(temporal_filter.vecteur_ratio(ll+1)*end):end,:),[],3)];
                    
                end
            end
        end
        
        % desallocation avant reallocation
        
        clear trajectory_subsets  kspace_subsets
        
        kspace_subsets=kspace_subsets_with_filter;
        trajectory_subsets=trajectory_subsets_with_filter;
        
    else        
        
    end
    
    
    for j=1:1:number_of_output_respiratory_phases
        for i=1:number_of_output_cardiac_phases
            
            t = cputime;
            
            %% DCF
            
            disp('starting density compensation');
            info.verbose = 0;
            info.numIter = 1;
            info.effMtx  = enc_Nx;
            DCF = sdc3_MAT(reshape(trajectory_subsets{i,j},[],3)',info.numIter,info.effMtx,info.verbose,info.osf);                       
            
            %% gpuNUFFT
            disp('starting gpuNUFFT');
            info.wg =2;
            info.sw = 4;
            mattest=reshape(trajectory_subsets{i,j},[],3)';
            FT = gpuNUFFT(mattest,DCF,info.osf,info.wg,info.sw,[info.N,info.N,info.nSl],[],true);
                                    
            %% Coil combine
            
            disp('starting coil combine');            
            image_nocrop(:,:,:) = coil_combine( kspace_subsets{i,j} , FT, info.N);    
            
            %% crop image and normalize
                       
            image_crop(:,:,:,i,j)  = crop_and_normalize_image( image_nocrop, 2 , info.N/2 , 2^12);
            
            e = cputime-t;  
            str_msg=sprintf(' ecg  %d  resp %d  timing  %f \n' ,i,j,e); disp(str_msg); 
            
        end
        
    end
    
end



%  filename_save = [ '/home/valery/Tempo/cs_traj.mat' ];
%  save(filename_save,'trajectory_subsets');
%  filename_save = [ '/home/valery/Tempo/cs_rawdata.mat' ];
%  save(filename_save,'kspace_subsets' , '-v7.3');

%
%     if (strcmp(mode.save_image,'Y'))
%
%         str_repetition_for_reco_indice=sprintf('%d',repetition_for_reco_indice-1);
%
%         filename_save = [ '/home/valery/Tempo/test_',str_repetition_for_reco_indice ,'.mat' ]
%
%         t = cputime;
%         disp('saving file');
%         save(filename_save,'image_crop');
%
%         e = cputime-t;
%         str_msg=sprintf('temps de sauvegarde %f \n',e ) ;
%
%     end

% plot3(trajectory_subsets(10,:,1),trajectory_subsets(10,:,2),trajectory_subsets(10,:,3),'*')

% close(figure(6))
% figure(6)
% for ii=1:10:size(image,3)
%     imagesc(image2(:,:,ii,1));  title( int2str(ii)); colormap(gray);
%     pause(0.5)
% end

% image_crop_save_data_ok=image_crop_save_data(:,:,:,:,1);
% image_crop_save_data=image_crop;
% filename_save = [ '/home/valery/Tempo/image_check.mat' ];
% save(filename_save,'image_crop_save_data_ok' );


% A=load(filename_save);
end



