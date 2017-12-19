
clear all
close all

folder_save='/home/valery/TempoFigures/';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=20;
% number_of_output_respiratory_phases=1;
% mode.filter='N';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=20;
% number_of_output_respiratory_phases=3;
% mode.filter='N';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=20;
% number_of_output_respiratory_phases=3;
% mode.filter='Y';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );


% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=1;
% number_of_output_respiratory_phases=3;
% mode.filter='N';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=1;
% number_of_output_respiratory_phases=5;
% mode.filter='N';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=1;
% number_of_output_respiratory_phases=7;
% mode.filter='N';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=1;
% number_of_output_respiratory_phases=10;
% mode.filter='N';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );
% 

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number_of_output_cardiac_phases=20;
% number_of_output_respiratory_phases=5;
% mode.filter='N';
% 
% [ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
% save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_of_output_cardiac_phases=20;
number_of_output_respiratory_phases=7;
mode.filter='N';

[ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_of_output_cardiac_phases=20;
number_of_output_respiratory_phases=10;
mode.filter='N';

[ image_crop ] = reco_UTE_dense_coil_ecg_filtre_function( number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter);
save_image_volume(image_crop, folder_save , number_of_output_respiratory_phases ,  number_of_output_cardiac_phases,  mode.filter );


 