%
% imOut=phaseOffsetCorrection()
%
% Author:   Aurélien TROTIER  (a.trotier@gmail.com)
% Date:     2016-06-22
% Partner: none
% Institute: CRMSB (Bordeaux)
%
% Function description:

%
% Input:
%       paramStructure : généré par la fonction readParams_Bruker.m
%       (potentiellement modifiée ensuite)
% Output:
%
%   phaseOff = exp(-i*2pi*phase)  à multiplier avec le kspace
%
%
% Algorithm & bibliography:
%
% See also :
%
% To do :
%   *fonctionne pour la 3D à modifier pour la 2D

function phaseOff = phaseOffsetCorrection(paramStructure)

sx=paramStructure.acqp.ACQ_size(1)/2;
sy=paramStructure.method.PVM_Matrix(2);


if(strcmp(paramStructure.method.PVM_SpatDimEnum,'<3D>'))
    sz=paramStructure.method.PVM_Matrix(3);
    phaseOffy=zeros(sx,sy,sz);
    phaseOffz=zeros(sx,sy,sz);
    
    for k=1:sz
        phaseOffz(:,:,k)=ones(sx,sy)*(k-sz/2+1)/sz*sz*paramStructure.method.PVM_SPackArrSliceOffset/paramStructure.method.PVM_Fov(3);
    end
else
    sz=1;
    phaseOffy=zeros(sx,sy,sz);
    phaseOffz=zeros(sx,sy,sz);
    phaseOffx=zeros(sx,sy,sz);
end




for k=1:sy
    phaseOffy(:,k,:)=ones(sx,sz)*(k-sy/2+1)/sy*sy*paramStructure.method.PVM_SPackArrPhase1Offset/paramStructure.method.PVM_Fov(2);
end

for k=1:sx
    phaseOffx(k,:,:)=ones(sy,sz)*(k-sx/2+1)/sx*sx*paramStructure.method.PVM_SPackArrReadOffset/paramStructure.method.PVM_Fov(1);
end

phaseOff=phaseOffz;
phaseOff=phaseOff+phaseOffx;
%phaseOff=phaseOff+phaseOffy;


phaseOff=exp(1i*2*pi*phaseOff);

end