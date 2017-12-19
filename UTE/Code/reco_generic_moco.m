clc;
close all;
clear all;

addpath('/home/valery/Reseau/Valery/MatlabGit/PCAMotionDescriptor_v02/');

%%=========================  Read data file ================================================


A=load('/home/valery/TempoFigures/image_resp_7_cardiac_1_filter_N.mat');
% 
volume=A.volume; clear A 

number_of_output_cardiac_phases=size(volume,4);




 for cardiac=1:number_of_output_cardiac_phases

volume_phase_1=double(squeeze(volume(41:168,21:148,1:128,cardiac,:)));

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

for resp = 1 : 3

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
  display_result2D(Iref(:,:,num_display_slice), ...
		   I(:,:,num_display_slice), ...
		   registered_image(:,:,num_display_slice), ...
		   motion_field(:,:,num_display_slice));

       
      registered_volume(:,:,:,cardiac,resp)=registered_image;
       
      
      
end


%%========================= Close the RealTItracker library ===========================

RTTrackerWrapper();

%%========================= Close the PCAMotionDescriptor library ===========================

PCAMotionDescriptorWrapper();


 end


 
 
volume_mean=mean(registered_volume,5);

view_cine( volume(41:168,41:168,41:168,:,1), volume_mean, number_of_output_cardiac_phases , 1, 1);



resp=2;
 for z=1:1:dimz
     
     subplot(1,3,1); imagesc(data_multislice(:, :, z, resp));
      subplot(1,3,2); imagesc(volume_mean(:, :, z, 1));
       subplot(1,3,3); imagesc(registered_volume(:, :, z, resp)-data_multislice(:, :, z, resp)); colormap(gray);
      pause(0.1);
 end


 
 figure(6)
  for jj=1:1:no_dyn
 subplot(2,3,jj); imagesc(squeeze(data_multislice(60, :, :, jj)));   colormap(gray);
  subplot(2,3,jj+no_dyn); imagesc(squeeze(registered_volume(60, :, :, jj)));   colormap(gray);
  end
  
  
  
  
 figure(6)
  for jj=1:1:no_dyn
 subplot(2,3,jj); imagesc(squeeze(data_multislice(60, :, :, jj)));   colormap(gray);
  subplot(2,3,jj+no_dyn); imagesc(squeeze(registered_volume(60, :, :, jj)));   colormap(gray);
  end
  
  
% %%========================= Performing the PCA ===========================
% disp('Compute fast PCA...');
% tStart = tic;
% 
% V0 = zeros(2*dimx*dimy,dimz,no_dyn_learning_step);
% 
% for z = 1 : dimz
%   [U0, S0, V0_] = svd(learned_motion(:,:,z),'econ');
%   V0(:,z,:) = V0_;
% end
% 
% tEnd = toc(tStart);
% fprintf('PCA done in %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
% 
% %%========================= Registration loop over the dynamically acquired images during the interventional procedure ======
% 
% %% Define registration parameters
% PCAMotionDescriptorWrapper(Iref, dimx, dimy, dimz, ...
% 			   id_registration_method, ...
% 			   V0,base_size, ...
% 			   nb_raffinement_level, ...
% 			   accelerationFactor, ...
% 			   alpha);
% 
% curve_to_plot = zeros(no_dyn,base_size);

% for i = 1 : no_dyn
%   
%   %% Get the current image  
%   I = data_multislice(:, :, :, i);
%   
%   %% Estimate the motion between the reference and the current images
%   [motion_descriptor, motion_field] = PCAMotionDescriptorWrapper(I);
%  
%   %% Apply the estimated motion on the current image
%   [registered_image] = RTTrackerWrapper(I,motion_field);
%   
%   %% Store the estimated motion descriptor
%   curve_to_plot(i,:) = motion_descriptor(:,num_display_slice);
%   
%   %% Display registered images & estimated motion field
%   display_result2D(Iref(:,:,num_display_slice), ...
% 		   I(:,:,num_display_slice), ...
% 		   registered_image(:,:,num_display_slice), ...
% 		   motion_field(:,:,num_display_slice));
% 
%   % Display motion descriptor  
%   figure(2);
%   plot(curve_to_plot); title('Estimated motion descriptor');
% 
% end


