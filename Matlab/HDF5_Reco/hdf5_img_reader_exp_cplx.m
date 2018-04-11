close all
clear all


[ str_user ] = get_PC_name();

 filename = ['/home/', str_user, '/Dicom/DIXON/Verification/Ex_Vivo/2D/20180326_Unmasc/31_33_img.h5'];
 hinfo    = hdf5info(filename);

 res_m    = single(h5read(filename,  hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name));
 res_p    = single(h5read(filename,  hinfo.GroupHierarchy.Groups(1).Groups(2).Datasets(2).Name));  
 res_p    = (res_p-2048)/2048*pi; 
 
 res_cplx = res_m.*exp(1i.*res_p);
 exvivo   = permute(res_cplx,[1 2 5 3 4]);
 
 clear filename hinfo res_m res_p res_cplx
 
 
 filename = ['/home/', str_user, '/Dicom/DIXON/Verification/In_Vitro/2D/20180322/11_13_img.h5'];
 hinfo    = hdf5info(filename);

 res_m    = single(h5read(filename,  hinfo.GroupHierarchy.Groups(1).Groups(1).Datasets(2).Name));
 res_p    = single(h5read(filename,  hinfo.GroupHierarchy.Groups(1).Groups(2).Datasets(2).Name));  
 res_p    = (res_p-2048)/2048*pi; 
 
 res_cplx = res_m.*exp(1i.*res_p);
 invitro   = permute(res_cplx,[1 2 5 3 4]);
 
 clear str_user filename hinfo res_m res_p res_cplx