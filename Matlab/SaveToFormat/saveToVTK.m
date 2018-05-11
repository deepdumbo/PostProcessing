function saveToVTK(array, filename)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    dispstat('','init');
    
    if(isempty(regexp(filename,'*.vtk','ONCE')))
        tmp = [filename '.vtk'];
        filename = tmp;
        clear tmp;
    end

    [nx, ny, nz] = size(array);
    
    dispstat(sprintf(' - Matrix to VTK (size : %dx%dx%d)',nx,ny,nz),'keepthis','timestamp');
    
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
   % fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
   % fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
    %fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, 'SCALARS ImageFile float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    %fprintf(fid, '\n');
    

    for c=1:nz
        for b=1:ny
            for a=1:nx            
                fprintf(fid, '%.3f\t', array(a,b,c));
            end;
        end;
         dispstat(sprintf(' - Progress : %d/%d (%d%%)', c, nz, floor(100 * c/nz)));
    end;
    
    fclose(fid);
    
    dispstat(sprintf(' - Done.'),'keepprev','timestamp');
return