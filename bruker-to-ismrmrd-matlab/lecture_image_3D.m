close all
clear all
% filename_image='/home/valery/Dev/Bruker_data/out.h5'
% h5disp(filename_image)
% info = h5info(filename_image)
% name = info.Groups(1).Name
% fullname= [name,];
% 
% data = h5read('/home/valery/Dev/Bruker_data/out.h5', '/lala');

% data = h5read('/home/valery/Dev/Bruker_data/new.h5', '/data');
% data = h5read('/home/valery/Dev/Bruker_data/raw_bruker_5.h5', '/dataset/data');


% clear data
% data = h5read('/home/valery/Dev/Bruker_data/out.h5','/Image1/data')

% addpath('/home/valery/Reseau/Valery/MatlabGit/Generic_functions/');


file='/home/valery/Dev/gadgetron-ihu/build/out_00000.real'

nX=200;
nY=200;
nZ=220;

A=read_binary_array(file,3, [nX ,nY ,nZ] );

close(figure(1))
figure(1)
for z = 1:nZ
    
    imagesc(A(:,:,z)); colormap(gray);
    pause(0.1)
end