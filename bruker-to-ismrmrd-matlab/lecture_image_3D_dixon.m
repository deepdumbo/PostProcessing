close all
clear all


[ str_user ] = get_PC_name();

% [ str_network_imagerie, str_network_perso ] = get_network_name( str_user );

filename=['/home/',str_user,'/Reseau/Imagerie/For_Valery/TestRecoKylian/ExVivo/FID/image_fatwater_39.h5'];
hinfo = hdf5info(filename);


res_m=single(h5read(filename,  hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name));
res_p=single(h5read(filename,  hinfo.GroupHierarchy.Groups(1).Groups(2).Datasets(2).Name));

res_p=(res_p-2048)/2048*pi;

% res_cplx=res_m.*exp(1i.*res_p);


figure(1)
for s=1:10
    subplot(2,5,s); imagesc(res_m(:,:,1,1,s)); colormap(gray);
end

figure(2)
for s=1:2
    subplot(2,3,s); imagesc(res_p(:,:,1,1,s)); colormap(gray);
end
