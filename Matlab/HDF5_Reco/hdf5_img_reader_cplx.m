function [data, hinfo] = hdf5_img_reader_cplx(data_loc)

% Get the info from .h5 file
%acq_path = '/home/kylianh/Dicom/DIXON';
%data_loc = fullfile(acq_path, img_name);
fprintf(1, '\nImage data directory :\n%s\n', data_loc);

fprintf(1, '\nRead hinfo from .h5 file.\n');
hinfo = hdf5info(data_loc);

% Extract the data

fprintf(1, '\nExtract data from :\n%s\n', hinfo.GroupHierarchy.Groups(1).Groups.Datasets(2).Name);
data_struct = h5read(data_loc, hinfo.GroupHierarchy.Groups(1).Groups.Datasets(2).Name);

% Create a complex based on the data_struct
fprintf(1, '\nCreate complex matrix from data structure.\n');
data(:,:,:) = complex(data_struct.real(:,:,:,1,:), data_struct.imag(:,:,:,1,:));

end