clc;
close all;
clear all;

addpath(pwd);

%%=========================  Read data file ================================================

number_of_output_respiratory_phases=3;
number_of_output_cardiac_phases=20;

for method=1:1
    
    for j=1:1:number_of_output_respiratory_phases
        for i=1:number_of_output_cardiac_phases
                        
            str_i=sprintf('%d',i );
            str_j=sprintf('%d',j );
            
            fprintf('method: %d resp: %d cardiac: %d\n',method, j, i );
            
            if (method==1)
                str_method='adapt';
            elseif (method==2)
                str_method='inati';
            else
                str_method='';
            end
            
            filename_img_cgsense=['/home/valery/Tempo/' , 'img_cgsense_', str_j , '_', str_i , '_',str_method , '.mat'];
            filename_img_comb=['/home/valery/Tempo/' , 'img_comb_', str_j , '_', str_i , '_',str_method ,  '.mat'];
            
            A=  load(filename_img_comb);
            B=  load(filename_img_cgsense);
            
            if (method==1)                
                img_comb_adapt(:,:,:,i,j)=abs(A.img_comb);
                img_cgsense_adapt (:,:,:,i,j)=abs(B.img_cgsense); 
            end  
            
        end
    end
end


% A=load('/home/valery/TempoFigures/image_resp_3_cardiac_20_filter_N.mat');

volume_2=img_cgsense_adapt; 

number_of_output_cardiac_phases=size(volume_2,4);

size(volume_2)

x_min=41;
x_max=168;
y_min=41;
y_max=168;
z_min=41;
z_max=168;

 for cardiac=1:number_of_output_cardiac_phases


volume_phase_1=double(squeeze(volume_2(x_min:x_max,y_min:y_max,z_min:z_max,cardiac,:)));

size(volume_phase_1)

% [data,dimx,dimy,dimz,no_dyn] = load_dat('/home/valery/Reseau/Valery/MatlabUnix/2015-06-08-OpticalFlow_for_Solenn_Valery/data/abdomen2D.dat');

no_dyn=size(volume_phase_1,4);
dimx=size(volume_phase_1,1);
dimy=size(volume_phase_1,2);
dimz=size(volume_phase_1,3);
%%========================= Configuration parameters for the motion estimation library =====

%% Define registration method
%% 0: No motion estimation
%% 1: L2L2 optical flow algorithm
%% 2: L2L1 optical flow algorithm
id_registration_method = 2;

% Dynamic image used as the reference position
reference_dynamic = 2; 

%% Weighting factor (between 0 and 1) Close to 0: motion is highly sensitive to grey level intensity variation. Close to 1: the estimated motion is very regular along the space. See http://bsenneville.free.fr/RealTITracker/ for more informations
alpha = 0.3;   
if (id_registration_method == 2)
  alpha = 0.6;
end

%% Computation of the highest resolution level to perform
%% (accelerationFactor=0 => computation is done on all resolution levels,
%%  accelerationFactor=1 => computation is done on all resolution levels except the highest one)
accelerationFactor = 0;

%% Number of iterative raffinement within each resolution level 
nb_raffinement_level = 3;    

%% Number of dynamic in the learning step
no_dyn_learning_step = 20;

%% Number of coefficient in the motion descriptor
base_size = 5;

%% Optional switch to distinguish 2D multislice registration/3D registration (if no parameter is set, a 2D registration is performed for dimz=1 and a 3D registration is performed for dimz>1)
do_2D_registration = 0;

%% Select the slice for which we will display the results
num_display_slice = 2;

%%========================= Adjustement of grey level intensities =========================


% %% Get the reference image for the registration
Iref = volume_phase_1(:, :, :, reference_dynamic);

% %% Normalize the reference image
aux = Iref;
Iref = (aux - min(aux(:)))/(max(aux(:)) - min(aux(:)));


%% Normalize all other images by adjusting the mean of the images (less sensitive to local grey level variations compared to a "min-max" method)
for resp = 1 : no_dyn
  
    aux_ref = Iref;
    aux = volume_phase_1(:, :, :, resp);
    data_multislice(:, :, :, resp) = aux * (mean(aux_ref(:))/mean(aux(:)));
  
end

%%========================= Initialisation of the RealTItracker library =============

%% Define registration parameters
RTTrackerWrapper(dimx, dimy, dimz, ...
		 id_registration_method, ...
		 nb_raffinement_level, ...
		 accelerationFactor, ...
		 alpha, ...
		 do_2D_registration);  

%%========================= Registration loop over the dynamically acquired images during the learning step ======

learned_motion = zeros(no_dyn_learning_step,2*dimx*dimy,dimz);

for resp = 1 : no_dyn

  %% Get the current image  
  I = data_multislice(:, :, :, resp);
  
  %% Estimate the motion between the reference and the current images
  RTTrackerWrapper(Iref, I);
  
  % Apply the estimated motion on the current image
  [registered_image] = RTTrackerWrapper(I);
  
  % Get the estimated motion field
  [motion_field] = RTTrackerWrapper();

%   learned_motion(i,1:dimx*dimy,:) = reshape(real(motion_field),dimx*dimy,dimz);
%   learned_motion(i,dimx*dimy+1:2*dimx*dimy,:) = reshape(imag(motion_field),dimx*dimy,dimz);

  %% Display registered images & estimated motion field
%   display_result2D(Iref(:,:,num_display_slice), ...
% 		   I(:,:,num_display_slice), ...
% 		   registered_image(:,:,num_display_slice), ...
% 		   motion_field(:,:,num_display_slice));
       
      registered_volume(:,:,:,cardiac,resp)=registered_image;       
      
      
end


%%========================= Close the RealTItracker library ===========================

RTTrackerWrapper();

%%========================= Close the PCAMotionDescriptor library ===========================

PCAMotionDescriptorWrapper();


 end


 
 
registered_volume_mean=mean(registered_volume,5);

volume_2_mean=mean(img_comb_adapt(x_min:x_max,y_min:y_max,z_min:z_max,:,:),5);

% view_cine( volume_2(41:168,41:168,41:168,:,1), volume_mean, number_of_output_cardiac_phases , 1, 1);

addpath('/home/valery/Dev/cute/UTE/Code/')




view_cine_small( volume_2_mean,registered_volume_mean,  number_of_output_cardiac_phases , 1, 1);





cardiac=1;

% resp=2;

 for z=1:1:dimz
     
     subplot(1,2,1); imagesc(volume_2_mean(:, :, z, cardiac));
      subplot(1,2,2); imagesc(registered_volume_mean(:, :, z, cardiac));
%        subplot(1,3,3); imagesc(registered_volume_mean(:, :, z, 1)-data_multislice(:, :, z, resp)); colormap(gray);
      pause(0.15);
 end


 indice=50;
 figure(6)
  for resp=1:1:no_dyn
 subplot(2,no_dyn,resp); imagesc(squeeze(data_multislice(indice, :, :, resp)));   colormap(gray); title('no correction');
  subplot(2,no_dyn,resp+no_dyn); imagesc(squeeze(registered_volume(indice, :, :,end, resp)));   colormap(gray);  title('correction');
  end
  
  
  
 figure(6)
  for resp=1:1:no_dyn
 subplot(2,3,resp); imagesc(squeeze(data_multislice(60, :, :, resp)));   colormap(gray);
  subplot(2,3,resp+no_dyn); imagesc(squeeze(registered_volume(60, :, :, resp)));   colormap(gray);
  end
  
  


