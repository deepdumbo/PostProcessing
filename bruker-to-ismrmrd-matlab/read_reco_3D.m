

% https://bfs.u-bordeaux.fr/telecharge.php?choix=files/7a5d6ea48bf9d5c1b1895ebfb7b3e41f/testdata_heart_example_2_IPAT2.h5
% https://bfs.u-bordeaux.fr/telecharge.php?choix=files/3f8ef07179dac563736c91d33652a080/testdata_heart_example_1_IPAT2.h5
% https://bfs.u-bordeaux.fr/telecharge.php?choix=files/1adb4f71d928744dd4e4248cd7ff7c6b/testdata_heart_example_random_2.h5


clear all

filename_grappa='/home/valery/DICOM/out1_grappa_00000.real';

ndims=3;
dims=[128 128 64];
data_grappa_1 = read_binary_array(filename_grappa, ndims, dims);


filename_spirit='/home/valery/DICOM/out1_spirit_00000.real';

ndims=3;
dims=[128 128 64];
data_spirit_1 = read_binary_array(filename_spirit, ndims, dims);



size(data_grappa_1)
size(data_spirit_1)


figure(1)
for z=1:1:dims(3)
    
    subplot(121); imagesc(data_grappa_1(:,:,z)); colormap(gray);
    subplot(122); imagesc(data_spirit_1(:,:,z)); colormap(gray);
    pause(0.1);
    
end

figure(2)
z=40;
subplot(121); imagesc(data_grappa_1(:,:,z)); colormap(gray);
    subplot(122); imagesc(data_spirit_1(:,:,z)); colormap(gray);
    
  
    
   
    
filename_grappa='/home/valery/DICOM/out2_grappa_00000.real';

ndims=3;
dims=[128 128 64];
data_grappa_2 = read_binary_array(filename_grappa, ndims, dims);


filename_spirit='/home/valery/DICOM/out2_spirit_00000.real';

ndims=3;
dims=[128 128 64];
data_spirit_2 = read_binary_array(filename_spirit, ndims, dims);



figure(1)
for z=1:1:dims(3)
    
    subplot(121); imagesc(data_grappa_2(:,:,z)); colormap(gray);
    subplot(122); imagesc(data_spirit_2(:,:,z)); colormap(gray);
    pause(0.1);
    
end  
    
    
    
    
filename_grappa='/home/valery/DICOM/outh_grappa_00000.real';

ndims=3;
dims=[200 200 64];
data_grappa_humain = read_binary_array(filename_grappa, ndims, dims);


filename_spirit='/home/valery/DICOM/outh_spirit_00000.real';

ndims=3;
dims=[200 200 64];
data_spirit_humain = read_binary_array(filename_spirit, ndims, dims);



figure(1)
for z=1:1:dims(3)
    
    subplot(121); imagesc(data_grappa_humain(:,:,z)); colormap(gray);
    subplot(122); imagesc(data_spirit_humain(:,:,z)); colormap(gray);
    pause(0.1);
    
end




 

filename_spirit='/home/valery/Dev/gadgetron-ihu/build/out_00000.real';

ndims=3;
dims=[200 200 64];
data_spirit_random = read_binary_array(filename_spirit, ndims, dims);

figure(1)
for z=1:1:dims(3)
    
    subplot(121); imagesc(data_spirit_random(:,:,z)); colormap(gray);
%     subplot(122); imagesc(data_spirit_humain(:,:,z)); colormap(gray);
    pause(0.1);
    
end
   

