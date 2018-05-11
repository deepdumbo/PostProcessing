%% Convert data to ISMRM's Fat Water Toolbox format
clear imDataParams im_tb

  im = tmp;

% Considering 2D
if(ndims(im) == 4)
    [nX, nY, nCoils, nEchoes] = size(im)
    nZ = 1;
elseif(ndims(im) == 5)
    [nX, nY, nZ, nCoils, nEchoes] = size(im)
end

if(size(TE,2) == 1)
    imDataParams.TE = TE'/ 1000
else
    imDataParams.TE = TE / 1000
end
imDataParams.FieldStrength = 9.4
imDataParams.PrecessionIsClockwise = 1

im_tb = zeros(nX,nY,nZ,nCoils,nEchoes);
im_tb(:,:,:,:,:) = im;
imDataParams.images = im_tb;
save('VerifAllData/20180412_Unmasc2D_2_fw.mat','imDataParams');