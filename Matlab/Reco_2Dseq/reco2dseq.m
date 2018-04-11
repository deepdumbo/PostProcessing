function [im,rawData]=reco2dseq(acquisition_path, pdata_file, RecoMode)

[paramStructure]=readParams_Bruker('dirPath',acquisition_path);

    %Read file
    fid=fopen([paramStructure.dirPath '/pdata/' num2str(pdata_file) '/2dseq'],'r','b');
    rawData=fread(fid,'bit16','l');
    fclose(fid);
    
    % If 2D
    if(size(paramStructure.PVM_Matrix, 2) < 3)
        paramStructure.PVM_Matrix(3) = 1;
    end
    
    % If coils combined or not
    if(strcmpi(RecoMode,'combined'))
        nCoils = 1;
    elseif(strcmpi(RecoMode,'shuffle'))
        nCoils = 7;
    end
    
im= reshape(rawData,paramStructure.PVM_Matrix(1),paramStructure.PVM_Matrix(2),paramStructure.PVM_Matrix(3),nCoils,[]);
    
end