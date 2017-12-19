function  data = read_binary_array(filename, ndims, dims)
%  data = read_gadgetron_array(filename)
%  
%  Reads simplified array format output from the Gadgetron
%
%  The datatype is determined by the file extension.
%     - *.short : 16-bit unsigned integer
%     - *.real  : 32-bit float
%     - *.cplx  : 32-bit complex (two 32-bit values per data element)
%
%
if (~exist(filename,'file')),
    error('File not found.');
end

[path,name,ext] = fileparts(filename);

ext = lower(ext);

if (~strcmp(ext,'.short') && ~strcmp(ext,'.real') && ~strcmp(ext,'.cplx') && ~strcmp(ext,'.bin')&& ~strcmp(ext,'.bind')),
   error('Unknown file extension'); 
end



f = fopen(filename);
%ndims = fread(f,1,'int32'); 
%dims = fread(f,ndims,'int32'); 

switch ext
    case '.short'
        data = fread(f,prod(dims),'uint16'); 
    case '.real'
        data = fread(f,prod(dims),'float32'); 
    case '.bin'
        data = fread(f,prod(dims),'float32');  
    case '.bind'
        data = fread(f,prod(dims),'double');     
    case '.cplx'
        data = fread(f,2*prod(dims),'float32'); 
        data = complex(data(1:2:end),data(2:2:end));
    otherwise     
end

fclose(f);

%original=data;
data = reshape(data,dims);

end



