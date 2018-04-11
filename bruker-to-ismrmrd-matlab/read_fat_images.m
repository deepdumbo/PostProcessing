

% Get the info from .h5 file
%acq_path = '/home/kylianh/Dicom/DIXON';
%data_loc = fullfile(acq_path, img_name);

data_loc='/home/valery/Dev/Data/out.h5'

fprintf(1, '\nImage data directory :\n%s\n', data_loc);

fprintf(1, '\nRead hinfo from .h5 file.\n');
hinfo = hdf5info(data_loc);

% Extract the data

fprintf(1, '\nExtract data from :\n%s\n',hinfo.GroupHierarchy.Groups.Groups(1).Datasets(2).Name);
res1 = h5read(data_loc, hinfo.GroupHierarchy.Groups.Groups(1).Datasets(2).Name);

fprintf(1, '\nExtract data from :\n%s\n', hinfo.GroupHierarchy.Groups.Groups(2).Datasets(2).Name);
res2 = h5read(data_loc, hinfo.GroupHierarchy.Groups.Groups(2).Datasets(2).Name);

figure()
for s=1:4
    subplot(2,2,s)
imagesc(res1(:,:,1,1,s)); colormap(gray);
end



figure()
for s=1:2
    subplot(2,2,s)
imagesc(res2(:,:,1,1,s)); colormap(gray);
end