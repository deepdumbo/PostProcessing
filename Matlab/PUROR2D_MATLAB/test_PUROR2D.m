%% % % % % % % % % % % % % % % % % % % % % 
%       Junmin Liu (junminliumri@gmail.com)
%       University of Western Ontario, Canada
%       Feb 04, 2014
% % % % % % % % % % % % % % % % % % % % %
% modified by Junmin liu Feb 04, 2014
% PUROR 3D and MEX files are also avalable upon request 
%
%   Reference:
%   ---------
%   Junmin Liu and Maria Drangova, 
%"Intervention-based multidimensional phase unwrapping using recursive orthogonal referring",
%   Magnetic Resonance in Medicine, Volume 68(4):1303-–1316, 2012
%
%
clear all;
close all;
load testPUROR.mat;
%
debug_PUROR = 0;
%----------------------------------------------------------------------
% determie the threshold fot noise 
[th4noise] = NoiseEstimation(complex_image);
% prepare the mask
th4unwrap = 0.0; % the pixels included in the unwrapping mask 
th4supp = 4.0; % the pixels included in the support mask
th4stack = 8.0; % the pixels used for stacking the multiple slices
[mask4unwrap, mask4supp, mask4stack] = mask_generation(complex_image,th4noise,th4unwrap,th4supp,th4stack);
%---------------------------------------------------------------------
[ unwrapped_phase ] = PUROR2D( complex_image,mask4unwrap, mask4supp, mask4stack, debug_PUROR );
%---------------------------------------------------------------------
matrix_size = size(unwrapped_phase);
if length(matrix_size) == 2
matrix_size(3) = 1;
end
%
close all
for index_slice = 1:matrix_size(3)
   figure(1)
   imagesc(unwrapped_phase(:,:,index_slice));colormap gray;axis square;axis off;
   index_slice
   pause
end

