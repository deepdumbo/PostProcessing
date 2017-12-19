function [ trajRAD ] = ComputeUTEtraj( N,nFE,SGPoints,osf,BWp ,RiseT ,GA)
%COMPTEUTETRAJ


sizeR2=N/2*osf;
trajBase=ones(sizeR2,1);

BWt=BWp*sizeR2*2;
Dt=1/BWt; %dwell time

NbPointRamp=RiseT*(10^-6)/Dt; %# points in ramp


for i=1:SGPoints
    trajBase(i)=0;
end

for i=0:floor(NbPointRamp)-1
   trajBase(i+1+SGPoints)=i/NbPointRamp;

end

a=cumtrapz(trajBase);
trajx = [];
trajy = [];
trajz = [];
trajRAD=[];

x=a;
phin1=0;

for ii=1:nFE
    
    hn=-1+2*(ii-1)/(nFE-1);
    theta=acos(hn);
    
    if ii==nFE || ii==1
        phi=0;
    else
        phi=mod(phin1+3.6/sqrt(nFE*(1-hn^2)),2*pi);
    end
    
    phin1=phi;
    
    xA = cos(phi)*sin(theta);
    yA = sin(phi)*sin(theta);
    zA = cos(theta);
    
    trajx(1:sizeR2,ii)=x*xA;
    trajy(1:sizeR2,ii)=x*yA;
    trajz(1:sizeR2,ii)=x*zA;
    
end

trajxRAD=trajx;
trajyRAD=trajy;
trajzRAD=trajz;

trajxRAD=reshape(trajxRAD,[],1);
trajyRAD=reshape(trajyRAD,[],1);
trajzRAD=reshape(trajzRAD,[],1);

trajRAD(:,1) = trajxRAD(:);
trajRAD(:,2) = trajyRAD(:);
trajRAD(:,3) = trajzRAD(:);

trajRAD=trajRAD/(2*max(max(abs(trajRAD))));
trajRAD=-trajRAD;

if GA==2
    trajRAD2=reshape(trajRAD,[N/2,nFE,3]);
    size(trajRAD2)
    temp=trajRAD2;
    nFE
    for i=1:nFE
        if i==nFE
            index=1;
        else
            index=ceil(mod((nFE-(i))*((i)*(sqrt(5)-1)/2),double(nFE-(i))));
        end
        
        trajRAD2(:,i,:)=temp(:,index,:);
        temp=permute(temp,[1,3,2]);
        temp(:,:,index)=[];
        temp=permute(temp,[1,3,2]);
        trajRAD=reshape(trajRAD2,[N/2*nFE,3]);
        i
    end
    
end

end

