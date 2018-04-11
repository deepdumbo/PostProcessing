close all
clear all

[ str_user ] = get_PC_name();

filename=['/home/',str_user,'/Dev/Data/out_cmrr_matlab.h5'];
hinfo = hdf5info(filename);

res.magnitude=h5read(filename, hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name);
res.phase=h5read(filename, hinfo.GroupHierarchy.Groups(1).Groups(2).Datasets(2).Name);

close(figure(1))
figure(1)

for s=1:size(res,3)
subplot(131); imagesc(res(:,:,s,:,1));  colormap(gray); 
subplot(132); imagesc(res(:,:,s,:,2)); title(num2str(s)); colormap(gray);
subplot(133); imagesc(res(:,:,s,:,1)-res(:,:,s,:,2));  colormap(gray);
pause(0.10)
end


close(figure(1))
figure(1)
s=60;
subplot(131); imagesc(res(:,:,s,:,1)); title('echo1'); colormap(gray);
subplot(132); imagesc(res(:,:,s,:,2)); title('echo2'); colormap(gray);
subplot(133); imagesc(res(:,:,s,:,1)-res(:,:,s,:,2)); title('substraction');  colormap(gray);
