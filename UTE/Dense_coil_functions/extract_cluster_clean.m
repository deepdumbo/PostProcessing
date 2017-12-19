function [ cluster_sum , coilID] = extract_cluster_clean( matricG, matrixC, input, thresholdU )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nbPointsZero=size(input,1);
nCh=size(input,2);

subgroup = zeros(nCh,1);
subgroup(abs(matricG(:,1))>thresholdU)=1;
subindex = find(subgroup==1);
coilID = subgroup;

count=0;
cluster_sum=zeros(size(input,1),1);

for c=1:nCh
    
    if (abs(matricG(c,1))>thresholdU)
        
        %                 msg_str=sprintf('npZ: %d j :%d c : %d   %f \n',npZ, j, c,UUmatrix(c,j,npZ) );  disp(msg_str);
        count=count+1;
        
        if (matrixC(c,1)>0)
            cluster_sum=cluster_sum+input(:,c);
        else
            cluster_sum=cluster_sum-input(:,c);
        end
        %                              subplot(2,1,j+(npZ-1)*2);  plot(t,squeeze(inputdata(npZ,c,:)));    hold on;  xlim([0, size(inputdata,3)*T]);
    end;
    
    %             subplot(2,2,j+(npZ-1)*2);plot(t',cluster_sum,'r'); xlim([0, size(inputdata,3)*T]);
end

if (count>1)
    cluster_sum=cluster_sum/count;
end

% msg_str=sprintf('count %d \n',count);
% title(msg_str);


if (count>1)

subindex'

end

end

