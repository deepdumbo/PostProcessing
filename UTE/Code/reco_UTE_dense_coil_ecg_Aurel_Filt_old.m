
%% DONE
% le script gere le self gating respiratoire
% le script gere le self gating cardiaque
% le script gere le self gating respiratoire et cardiaque

%% TODO

% la detection de peak cardiaque peut ??tre am??lior??e
% reechantionner l'ecg le vrai
% ajouter la MOCO


clear all
close all

addpath('../Generic_functions/');

[ str_user ] = get_PC_name();

folder= ['/home/', str_user , '/Reseau/Imagerie/DICOM_DATA/2017-03-08-UTE/';]
name = 'meas_MID00393_FID16077_UTE_SG1_5_postinf_FA30_NR10_2.h5';

filename = [ folder, '/FID/'  name ] ;

if (strcmp(name,'meas_MID00393_FID16077_UTE_SG1_5_postinf_FA30_NR10_2.h5'))
    folder='/home/valery/';
    filename = [ folder,  name ] ;
end

if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

folder.data='../../Data/';




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




%% quelques v??rifications

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

nbProjections=nbProjections;
readout=enc_Nx;

nbChannels=30;
channels=8;
repetition=10;


%% parametres du self gating

threshold.exclusion=1000;
threshold.debut=threshold.exclusion+1;
threshold.fin=nbProjections-threshold.exclusion;

threshold.resp=1;
threshold.ecg=1;
threshold.svd=0.1;

threshold.exclusion_all=1000;
threshold.debut_all=threshold.exclusion_all+1;
threshold.fin_all=nbProjections*repetition-threshold.exclusion_all;

threshold.exclusion_ecg=1000;
threshold.debut_ecg=threshold.exclusion_ecg+1;
threshold.fin_ecg=nbProjections*repetition-threshold.exclusion_ecg/2;


%% chargement des donn??es de self gating

data.input.raw=zeros(nbProjections,nbChannels,repetition);

for r=1:1:repetition
    
    str_r=sprintf('%d',r-1);
    
    filename=[folder.data,'/Navigator/dense_data_store_abs_',str_r,'.dat']
    temp=load(filename);
    data.input.raw(:,:,r)=temp;
    
end


%% on tri les donn??es et on applique le filtre sur toutes les projections

data.allrep.raw=reformat_to_one_block(data.input.raw);

data.allrep.filter.resp=apply_butterworth_filtering(data.allrep.raw,1, Fs_ute);

data.allrep.substract.resp=zeros(size(data.allrep.raw));

for c=1:1:nbChannels
    data.allrep.substract.resp(:,c)=data.allrep.raw(:,c)-data.allrep.filter.resp(:,c);
end

%% on reformate les donn??es par r??p??tition

[ data.filter.resp ] = reformat_to_each_repitition( data.allrep.filter.resp, repetition );

[ data.substract.resp ] = reformat_to_each_repitition( data.allrep.substract.resp, repetition );

%% ici on boucle sur les repetitions

[ navigator ] = extract_navigator_resp_and_ecg( data , threshold.debut, threshold.fin ,repetition , threshold , Fs_ute );

navigator.allrep.resp.normalize=navigator.resp.normalize(:);

navigator.allrep.ecg.filter=navigator.ecg.filter(:);


%% partie respiration

nb_phase_respiratoire=3;

[ table ] = extract_liste_respiratory( navigator.allrep.resp.normalize, threshold.debut_all, threshold.fin_all, nb_phase_respiratoire );


%% partie ecg

%% on prend manuellement les coils 21 et 22 car sur les repetitions 8 et 9 l'auto dectection ne marche pas

navigator.ecg.test=  (data.allrep.substract.resp(:,21)+ data.allrep.substract.resp(:,22))/2;

[ navigator.ecg.testfilter(:,1) ]=apply_butterworth_filtering(navigator.ecg.test(:,1),3, Fs_ute);



%% detection des peaks

[ data_peak_dectection.plus, data_peak_dectection.moins ] = find_peaks(navigator.ecg.testfilter, threshold.debut_ecg , threshold.fin_ecg);

[ time_difference.plus ] = find_time_difference_between_peaks( data_peak_dectection.plus, info.Tr_ );

[ time_difference.moins ] = find_time_difference_between_peaks( data_peak_dectection.moins, info.Tr_ );

% ca c'est just pour enlever le delai entre les pics + et -
delay=data_peak_dectection.plus(:,1)-data_peak_dectection.moins(:,1);

data_peak_dectection.plusdelay(:,1)=data_peak_dectection.plus(:,1)-median(delay);
data_peak_dectection.plusdelay(:,2)=data_peak_dectection.plus(:,2);

%% on enleve les peaks qui ont mal ??t?? d??tect??


threshold.repolarisation_time=300;
threshold.maximum_time=800;


data_peak_dectection.plus_new=remove_bad_peaks(data_peak_dectection.plus , time_difference.plus , threshold.repolarisation_time, info.Tr_ , threshold.maximum_time);

[ time_difference.plus_new ] = find_time_difference_between_peaks( data_peak_dectection.plus_new, info.Tr_ );

data_peak_dectection.moins_new=remove_bad_peaks(data_peak_dectection.moins , time_difference.moins , threshold.repolarisation_time, info.Tr_ , threshold.maximum_time);

[ time_difference.moins_new ] = find_time_difference_between_peaks( data_peak_dectection.moins_new, info.Tr_ );




repetition_for_reco=10;

%% sort ecg data

nb_phase_ecg=60 %g.NumberofPhase    %NB of phases

[ table_ecg ] = extract_list_ecg( data_peak_dectection.moins_new(:,1), nb_phase_ecg ,repetition_for_reco, nbProjections  );

debut_change=5;

[ table_ecg_in_systole ] = start_list_ecg_in_systole( table_ecg , debut_change, nb_phase_ecg);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lecture du vrai ecg

clear ecg_read ecg_truth ecg_tempo
% ecg_truth=load('/home/valery/Reseau/Imagerie/DICOM_DATA/2017-03-08-UTE/ecg/UTEPMUsignal_20170309T134113.mat');

filename=['/home/', str_user , '/Reseau/Imagerie/DICOM_DATA/2017-03-08-UTE/ecg/UTEPMUsignal_20170309T',info.study_time_read,'.mat'];

[ ecg_truth ] = read_ecg_from_mat( filename );


T_ecg = 1/400;            % Sample time
L_ecg = size(ecg_truth,1);  % Length of signal
t_ecg = (0:L_ecg-1)*T_ecg;        % Time vector

%t_ecg(end)/300000*1000

T_ute = t_ecg(end)/(nbProjections*10)   ; %1167/300000*1000; % Sample time
L_ute = size(navigator.ecg.testfilter,1);  % Length of signal
t_ute = (0:L_ute-1)*T_ute;        % Time vector


peak_detection.plus_true_timing(:,1)=data_peak_dectection.plus(:,1)*T_ute;
detection.moins_true_timing(:,1)=data_peak_dectection.moins(:,1)*T_ute;

peak_detection.plus_true_timing(:,2)=data_peak_dectection.plus(:,2);
peak_detection.moins_true_timing(:,2)=data_peak_dectection.moins(:,2);


[ ecg_resample ] = resample_ecg_to_ute_frequency( ecg_truth , nbProjections*10);



[ ecg_resample ] = resample_ecg_to_ute_frequency( ecg_truth , nbProjections*10);


[pks,locs]=findpeaks(-ecg_resample(:,2),'MinPeakDistance',10,'MinPeakHeight',200 );



close(figure(7))
figure(7)
subplot(1,1,1);
plot(t_ute,-ecg_resample(:,2),'r'); xlim([10 45]);
hold on;
plot(locs*T_ute, 200, '*');

[ table_ecg_recording ] = extract_list_ecg( locs(:,1), nb_phase_ecg ,repetition_for_reco, nbProjections  );

[ table_ecg_recording_extended ] = extract_list_ecg( locs(:,1), nb_phase_ecg*3 ,repetition_for_reco, nbProjections  );


x = diff(sign(ecg_resample(end:-1:1,2)));
indx_up = find(x>0);
test = find(x<0);

close(figure(7))
figure(7)
subplot(1,1,1);
plot(t_ute,ecg_resample(:,2),'r'); xlim([1034 1045]);
hold on;
plot(test*T_ute, 200, '*');


close(figure(7))
figure(7)
subplot(1,1,1);
plot(t_ecg,-ecg_truth(:,2)/2-400,'r');
hold on;
plot(t_ecg,ecg_truth(:,3)/2+200,'r');  xlim([1034 1045])
hold on;
plot(t_ecg,-ecg_truth(:,4)/10,'r');  %xlim([345 355])
hold on;
plot(t_ute,navigator.ecg.testfilter*2e5-600,'b');
hold on;
% plot(true_timing_peak_detection_plus(:,1),data_peak_dectection_plus(:,2), '*b');
hold on;
% plot(true_timing_peak_detection_moins(:,1),data_peak_dectection_moins(:,2), '*b');
hold on;
plot(data_peak_dectection.plus_new(:,1)*T_ute,data_peak_dectection.plus_new(:,2), '*b');
hold on;
plot(data_peak_dectection.moins_new(:,1)*T_ute,data_peak_dectection.moins_new(:,2), '*g');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





 





%%





clear  datatosave
datatosave(:,1)=t_ecg(:,1);
datatosave(:,2)=ecg_truth(:,3)/50;

min(datatosave(:,2))
max(datatosave(:,2))

save('/home/valery/Reseau/Valery/MatlabUnix/GPU/MakeFigures/ecg_voie2.txt','datatosave','-ascii');

datatosave(:,1)=t_ecg(:,1);
datatosave(:,2)=ecg_truth(:,4)/700;

min(datatosave(:,2))
max(datatosave(:,2))

save('/home/valery/Reseau/Valery/MatlabUnix/GPU/MakeFigures/ecg_trigger.txt','datatosave','-ascii');


clear datatosave
datatosave(:,1)=t_ute;
datatosave(:,2)=navigator.ecg.testfilter*2e4;

min(datatosave(:,2))
max(datatosave(:,2))

save('/home/valery/Reseau/Valery/MatlabUnix/GPU/MakeFigures/ecg_self.txt','datatosave','-ascii');


clear datatosave
datatosave(:,1)=data_peak_dectection.plus_new(:,1)*T_ute;
datatosave(:,2)=data_peak_dectection.plus_new(:,2);

min(datatosave(:,2))
max(datatosave(:,2))

save('/home/valery/Reseau/Valery/MatlabUnix/GPU/MakeFigures/ecg_self_trigger_plus.txt','datatosave','-ascii');


clear datatosave
datatosave(:,1)=data_peak_dectection.moins_new(:,1)*T_ute;
datatosave(:,2)=data_peak_dectection.moins_new(:,2);

min(datatosave(:,2))
max(datatosave(:,2))

save('/home/valery/Reseau/Valery/MatlabUnix/GPU/MakeFigures/ecg_self_triggermoins.txt','datatosave','-ascii');


%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% figure(4)
% plot(navigator.resp.normalize(debut_all:fin_all)); title('respiratory navigator'); xlim([1 30000]);
%
% close(figure(4))
% figure(4)
% plot(navigator.ecg.testfilter(debut_all:fin_all)); title('ecg navigator'); xlim([1 300000]);
% close(figure(6))
% figure(6)
% subplot(211)
% plot(table.respiration.affichage(:,1),'bo');
% hold on;
% plot(table.respiration.affichage(:,2),'co');
% hold on;
% plot(table.respiration.affichage(:,3),'go');
% hold on;
% plot(table.respiration.affichage(:,4),'yo'); ylim([0 1.5]);
% hold on;
% plot(table.respiration.affichage(:,5),'ro');
% hold on;
% plot(navigator.allrep.resp.normalize(debut_all:fin_all,1)); title('respiratory navigator'); xlim([1 30000]);
% subplot(212)
% plot(table.respiration_quartile.affichage(:,1),'bo');
% hold on;
% plot(table.respiration_quartile.affichage(:,2),'co');
% hold on;
% plot(table.respiration_quartile.affichage(:,3),'go');
% hold on;
% plot(table.respiration_quartile.affichage(:,4),'yo'); ylim([0 1.5]);
% hold on;
% plot(table.respiration_quartile.affichage(:,5),'ro');
% hold on;
% plot(navigator.allrep.resp.normalize(debut_all:fin_all,1)); title('respiratory navigator'); xlim([1 30000]);














close(figure(7))
figure(7)
subplot(1,1,1);
plot(t_ecg,-ecg_truth(:,2)/2-400,'r');
hold on;
% plot(t_ecg,ecg_truth(:,3)/2+200,'r');  xlim([345 355])
hold on;
% plot(t_ecg,-ecg_truth(:,4)/10,'r');
% hold on;
plot(t_ute,navigator.ecg.testfilter*2e5-600,'b');  xlim([100 120])
hold on;
% plot(true_timing_peak_detection_plus(:,1),data_peak_dectection_plus(:,2), '*b');
hold on;
% plot(true_timing_peak_detection_moins(:,1),data_peak_dectection_moins(:,2), '*b');
% hold on;
% plot(data_peak_dectection.plus_new(:,1)*T_ute,data_peak_dectection.plus_new(:,2), '*b');
hold on;
plot(data_peak_dectection.moins_new(:,1)*T_ute,data_peak_dectection.moins_new(:,2)-300, '*g');
hold on;
plot(t_ute, table_ecg*5-200);
hold on;
plot(t_ute, table_ecg_in_systole*5-300);

%% ici on pourrait supprimer tout ce qui n'est plus necessaire


%% lecture du kspace apres coil compression

% TODO chnager le nom de data input
% TODO afficher la liste avec le nombre de projection par repetition
repetition_for_reco=10;

input.matrice=zeros(readout, nbProjections ,channels, repetition_for_reco);

for r=1:1:repetition_for_reco
    
    t = cputime;
    
    str_r=sprintf('%d',r-1);
    
    filename= ['/home/', str_user ,'/Tempo/savekspace_ute_',str_r,'.bin']
    [ tempo ] =single( read_binary_complex( filename, [ readout*channels* nbProjections ] ));
    
    tempo1=reshape(tempo, [readout, channels,  nbProjections]);
    
    input.matrice(:,:,:,r)= permute(tempo1, [1 3 2]);
    
    e = cputime-t;
    str_msg=sprintf('%d %f  %f %f \n', r, e , (e -mod(e,60))/60 , mod(e,60)   ); disp(str_msg);
    
end

clear tempo tempo1

disp('fin lecture des donn??es');




% %% mode ecg
% mode.ecg_and_resp='Y';
% mode.ecg='N';
% mode.resp='N';
% mode.save_image='N';
%
% if (strcmp(mode.ecg_and_resp,'Y'));
%
%
%
%     size(table_input_resp)
%
%     table_input_ready=table_ecg_ok;
%
%     q=find(table_input_ready(:)>0);
%
%     size(q)
%
%     q=find(table_input_resp(:)~=1);
%
%     table_input_ready(q)=0;
%
%     q=find(table_input_ready(:)>0);
%
%     size(q)
%
%     nb_phase_en_tout=nb_phase_ecg;
%
% elseif (strcmp(mode.ecg,'Y'));
%
%     table_input_ready=table_ecg_ok;
%
%     size(table_input_ready)
%
%     nb_phase_en_tout=nb_phase_ecg;
%
% elseif (strcmp(mode.resp,'Y'));
%
%     table_input_resp=table.respiration.liste_ready(1:(nbProjections*repetition_for_reco_indice));
%
%     table_input_ready=table_input_resp;
%
%     disp('mode.resp');
%
% end


%% COMPUTE UTE traj


info.N=2*size(input.matrice,1);
info.nSl=info.N;

info.sample_time_us = 5; %% TODO a verifier et ?? rendre g??n??rique
info.osf       = 1;  %oversampling
info.SGPoints  = 2*acqOffSet;
info.GA        = AcqMode;
info.RiseT     = 300;
info.BWp       = (1/(info.sample_time_us*1e-6))/(info.osf*info.N); % Hz/px



% image_nocrop=zeros(info.N,info.N,info.N);
% image_crop=zeros(info.N/2,info.N/2,info.N/2,nb_phase_ecg, nb_phase_respiratoire);

clear  matrice_nRep  traj_nRep  rawDataECG trajECG FT mattest DCF matrice_nRep traj_nRep image_tempo_1 temp ktrajIdxAvant ktrajIdxApres ktrajIdx image_nocrop image_crop
clear rawDataECGtemp trajECGtemp

for repetition_for_reco_indice=10:1:repetition_for_reco;
    
    str_msg=sprintf(' reconstruction  %d / %d   \n', repetition_for_reco_indice ,  repetition_for_reco ); disp(str_msg);
    
    [enc_Nx,nbProjections,nCh,nRep_use]=size(input.matrice(:,:,:,1:repetition_for_reco_indice));
    
    str_msg=sprintf(' nPE  %d nCh %d  nFE %d  nRep_use %d \n', enc_Nx,nCh,nbProjections, nRep_use ); disp(str_msg);
    
    clear traj
    
    traj = ComputeUTEtraj(info.N, nbProjections , info.SGPoints, info.osf,info.BWp,info.RiseT, info.GA);
    
    clear matrice_nRep  traj_nRep
    
    matrice_nRep=input.matrice(:,:,:,1:repetition_for_reco_indice);
    
    %% Gather repetition
    % data
    matrice_nRep=permute(matrice_nRep,[1 3 2 4]);
    
    matrice_nRep=matrice_nRep(:,:,:);
    
    matrice_nRep=permute(matrice_nRep,[1 3 2]);
    
    str_msg=sprintf(' nPE  %d nCh %d  nFE %d  nRep_use %d \n', enc_Nx,nCh,nbProjections,nRep_use ); disp(str_msg);
    
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
    
    rawDataECG=cell(nb_phase_ecg,nb_phase_respiratoire);
    trajECG=cell(nb_phase_ecg,nb_phase_respiratoire);
    
    %% Image Reconstruction
    
    for j=1:1
        
        for i=1:nb_phase_ecg
            
            clear ktrajIdx
            
            %% on recupere la table de respiration
            [ table_input_ecg ] = remove_kspace_from_other_respiratory_stage(table.respiration.liste_ready(1:(nbProjections*repetition_for_reco_indice)) , j,  table_ecg_recording );
         
            
            %             size(q)
            
            ktrajIdx=find(table_input_ecg(:)==i);
            str_msg=sprintf(' phase resp %d  phase ecg  %d ,  nb projections %d ', j, i, size(ktrajIdx,1)); disp(str_msg);
            trajECG{i,j}=traj_nRep(:,ktrajIdx,:);
            rawDataECG{i,j}=matrice_nRep(:,ktrajIdx,:);
            
        end
        
    end
    
    
    %%  Realloc filt
    %je pars du principe qu'on a augment? le nombre de phaseECG par x3
    rawDataECGtemp=cell(length(2:3:nb_phase_ecg-1),nb_phase_respiratoire);
    trajECGtemp=cell(length(2:3:nb_phase_ecg-1),nb_phase_respiratoire);
    
    for j=1:1
        count=0;
        for i=2:3:nb_phase_ecg-1
            % Ajout des 3 phases les unes ? c?t? des autres
%             rawDataECGtemp{i,j}=[reshape(rawDataECG{i-1,j},[],1);reshape(rawDataECG{i,j},[],1);reshape(rawDataECG{i+1,j},[],1)];
%             trajECGtemp{i,j}=[reshape(trajECG{i-1,j},[],3);reshape(trajECG{i,j},[],3);reshape(trajECG{i+1,j},[],3)];
            count=count+1;
            
            fprintf('%d Complet : %d + %d + %d et 1/4 des hautes frequences de %d et %d\n',count, i-1,i,i+1,mod((i-1)-2,60)+1,mod((i-1)+2,60)+1);
            
            rawDataECGtemp{count,j}=[reshape(rawDataECG{i-1,j},[],1);reshape(rawDataECG{i,j},[],1)];
            trajECGtemp{count,j}=[reshape(trajECG{i-1,j},[],3);reshape(trajECG{i,j},[],3)];
            
            %Si on veut r?cup?rer une partie des donn?es des voisins.
            
            rawDataECGtemp{count,j}=[rawDataECGtemp{count,j};reshape(rawDataECG{mod((i-1)-2,60)+1,j}(round(3/4*end):end,:),[],1);reshape(rawDataECG{mod((i-1)+2,60)+1,j}(round(3/4*end):end,:),[],1)];
            trajECGtemp{count,j}=[trajECGtemp{count,j};reshape(trajECG{mod((i-1)-2,60)+1,j}(round(3/4*end):end,:,:),[],3);reshape(rawDataECG{mod((i-1)+2,60)+1,j}(round(3/4*end):end,:,:),[],3)];

        end
    end
    
    %%
    
    
    % desallocation
    
    clear traj_nRep matrice_nRep trajECG rawDataECG
    
    for j=1:1
        
        for i=1:size(rawDataECGtemp,1)
            
            t = cputime;
            
            %% DCF
            info.verbose = 0;
            info.numIter = 1;
            info.effMtx  = enc_Nx;
            DCF = sdc3_MAT(trajECGtemp{i,j}',info.numIter,info.effMtx,info.verbose,info.osf);
            
            disp('starting gpuNUFFT');
            
            %% gpuNUFFT
            info.wg =2;
            info.sw = 4;
            mattest=trajECGtemp{i,j}';
            FT = gpuNUFFT(mattest,DCF,info.osf,info.wg,info.sw,[info.N,info.N,info.nSl],[],true);
            
            disp('ending gpuNUFFT');
            
            %% Coil combine
            
            image_nocrop(:,:,:) = coil_combine( reshape (rawDataECGtemp{i,j}, []) , FT, info.N);
            
            e = cputime-t;
            
            str_msg=sprintf(' ecg  %d  resp %d  timing  %f \n' ,i,j,e); disp(str_msg);
            
            %% crop image
            osfcrop=2;
            Ncrop=info.N/2;
            image_tempo_1=image_nocrop(round((round(osfcrop*Ncrop)-Ncrop)/2+1):round(Ncrop+(round(osfcrop*Ncrop)-Ncrop)/2),round((round(osfcrop*Ncrop)-Ncrop)/2+1):round(Ncrop+(round(osfcrop*Ncrop)-Ncrop)/2),round((round(osfcrop*Ncrop)-Ncrop)/2+1):round(Ncrop+(round(osfcrop*Ncrop)-Ncrop)/2));
            
            %% IMAGE OUPUT
            image_crop(:,:,:,i,j)=(abs(image_tempo_1)-min(min(min(abs(image_tempo_1)))))* 2^12/(max(max(max(abs(image_tempo_1))))-min(min(min(abs(image_tempo_1)))));
            
        end
        
    end
    
    
end

%  filename_save = [ '/home/valery/Tempo/cs_traj.mat' ];
%  save(filename_save,'trajECG');
%  filename_save = [ '/home/valery/Tempo/cs_rawdata.mat' ];
%  save(filename_save,'rawDataECG' , '-v7.3');
 



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



% plot3(trajECG(10,:,1),trajECG(10,:,2),trajECG(10,:,3),'*')



% close(figure(6))
% figure(6)
% for ii=1:10:size(image,3)
%     imagesc(image2(:,:,ii,1));  title( int2str(ii)); colormap(gray);
%     pause(0.5)
% end


%image_crop_resp=image_crop;



close(figure(7))
figure(7)

for ii=50:2:size(image_crop,2)-50
    for jj=1:1:nb_phase_respiratoire
        subplot(2,3,jj)
        temp=squeeze(image_crop(:,ii,:,jj));
        imagesc(imrotate(temp,90));  title( int2str(ii)); colormap(gray);
    end
    pause(0.5)
end



close(figure(7))
figure(7)

for ii=124:2:124
    for jj=1:1:nb_phase_respiratoire
        subplot(2,3,jj)
        temp=squeeze(image_crop(ii,:,:,jj));
        imagesc(imrotate(temp,90));  title( int2str(ii)); colormap(gray);
    end
    pause(0.5)
end


close(figure(1))
figure(7)

for ll=74:2:74
    for jj=1:1:nb_phase_respiratoire
        for ii=1:1:nb_phase_ecg
            subplot(nb_phase_respiratoire,nb_phase_ecg,ii+ (nb_phase_ecg)*(jj-1));
            
            str_msg=sprintf('%d %d %d %d \n',ii,jj,ii+ (nb_phase_ecg-1)*jj ,jj); disp(str_msg);
            temp=squeeze(image_crop(ll,:,:,ii,jj));
            imagesc(imrotate(temp,90));  title( int2str(ii)); colormap(gray);axis square
        end
        
    end
    pause(0.5)
end







%% impact resp en coro

close(figure(8))
figure(8)
ii=2
indice=30;

while(1)
   for jj=1:1:nb_phase_respiratoire
        for kk=indice:10:130                    
                
                subplot(nb_phase_respiratoire,4, kk/10-2);            
                temp=squeeze(image_crop(:,kk,:,ii,jj));
                imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);        
        end
        
          pause(0.4)
   end        
   
end

%% impact resp en transverse

close(figure(9))
figure(8)
ii=2
indice=30;

while(1)
   for jj=1:1:nb_phase_respiratoire
        for kk=indice:10:130                    
                
                subplot(nb_phase_respiratoire,4, kk/10-2);            
                temp=squeeze(image_crop(:,:,kk,ii,jj));
                imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);        
        end
        
          pause(0.4)
   end        
   
end


%% impact resp en sagital

close(figure(9))
figure(8)
ii=2
indice=30;

while(1)
   for jj=1:1:nb_phase_respiratoire
        for kk=indice:10:130                    
                
                subplot(nb_phase_respiratoire,4, kk/10-2);            
                temp=squeeze(image_crop(kk,:,:,ii,jj));
                imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);        
        end
        
          pause(0.4)
   end        
   
end


close(figure(8))
figure(8)

taille=3;

indice=50;
while(1)
   
        for ii=1:1:16
             for jj=1:1:1
                 
%               for kk=indice:10:indice+40

%                 kk=50;
%                 subplot(nb_phase_respiratoire,4, kk-49 +(jj-1)*4);
%                 temp=squeeze(image_crop(:,kk,:,ii,jj));
%                 imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
                
                kk=70;
                subplot(nb_phase_respiratoire,taille, kk-69 +(jj-1)*taille);
                temp=squeeze(image_crop(kk,:,:,ii,jj));
                imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
                
                kk=80;
                subplot(nb_phase_respiratoire,taille, kk-79+1 +(jj-1)*taille);
                temp=squeeze(image_crop(:,kk,:,ii,jj));
                imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
                
                kk=90;
                subplot(nb_phase_respiratoire,taille, kk-89+2 +(jj-1)*taille);
                temp=squeeze(image_crop(:,kk,:,ii,jj));
                imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
%             end

              end
         pause(0.03)
        end          
end



close(figure(9))
figure(9)

for ii=20:1:size(image_crop,2)-20
    
    subplot(2,2,1)
    temp=squeeze(image_crop(:,ii,:,1));
    imagesc(imrotate(temp,90));  title( int2str(ii)); colormap(gray);
    
    subplot(2,2,2)
    temp=squeeze(image_crop(ii,:,:,1));
    imagesc(imrotate(temp,90));  title( int2str(ii)); colormap(gray);
    
    subplot(2,2,3)
    temp=squeeze(image_crop(:,:,ii,1));
    imagesc(temp);  title( int2str(ii)); colormap(gray);
    
    pause(0.3)
end





