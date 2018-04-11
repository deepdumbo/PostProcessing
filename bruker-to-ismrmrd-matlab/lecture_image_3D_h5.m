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






[ str_user ] = get_PC_name();

% [ str_network_imagerie, str_network_perso ] = get_network_name( str_user );

filename=['/home/',str_user,'/Dev/Data/out_spirit.h5'];
hinfo = hdf5info(filename);

for i=1:3
    
    res=h5read(filename, hinfo.GroupHierarchy.Groups(i).Groups.Datasets(2).Name);
    
    size(res)
    
    if(i==1)
        
        data_grappa_1=res;
        
    elseif (i==2)
        
        data_grappa_2=res;
        
    elseif (i==3)
                
        data_spirit_2=res;
        
    end
    
end

 [dimx, dimy, dimz, lala , echo ]= size(data_grappa_1);


figure(1)
for z=1:1:dimz
    
    subplot(121); imagesc(data_grappa_1(:,:,z)); colormap(gray);
    subplot(122); imagesc(data_grappa_2(:,:,z)); colormap(gray);
    pause(0.1);
    
end



% 
% 
% 
% % clear data
% % data = h5read('/home/valery/Dev/Bruker_data/out.h5','/Image1/data')
% % addpath('/home/valery/Reseau/Valery/MatlabGit/Generic_functions/');
% 
% file='/home/valery/Dev/gadgetron-ihu/build/out_00000.real'
% 
% nX=200;
% nY=200;
% nZ=220;
% 
% A=read_binary_array(file,3, [nX ,nY ,nZ] );
% 
% close(figure(1))
% figure(1)
% for z = 1:nZ
%     
%     imagesc(A(:,:,z)); colormap(gray);
%     pause(0.1)
% end