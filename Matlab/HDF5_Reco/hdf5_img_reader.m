function [data, hinfo] = hdf5_img_reader(data_loc)

% Get the info from .h5 file
fprintf(1, '\nImage data directory :\n%s\n', data_loc);

fprintf(1, '\nRead hinfo from .h5 file.\n');
hinfo = hdf5info(data_loc);

% Extract the data

fprintf(1, '\nExtract Magnitude data from :\n%s\n', hinfo.GroupHierarchy.Groups.Groups(1).Datasets(2).Name);
data_mag = single(h5read(data_loc, hinfo.GroupHierarchy.Groups.Groups(1).Datasets(2).Name));

fprintf(1, '\nExtract Phase data from :\n%s\n', hinfo.GroupHierarchy.Groups.Groups(2).Datasets(2).Name);
data_phase = single(h5read(data_loc, hinfo.GroupHierarchy.Groups.Groups(2).Datasets(2).Name));
data_phase = (data_phase-2048)/2048*pi;

% Create a complex based on the data_struct
fprintf(1, '\nCreate complex matrix from data structure.\n');
data(:,:,:) = data_mag(:,:,:,1,:).*exp(1i.*data_phase(:,:,:,1,:));

end