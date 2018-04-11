function [RSOS, noise, PHASE]=RSOS_function(data,no_rescale, noise)
     if nargin<2; no_rescale=0; end
%
% Convert k-space to RSOS
%

for n=1:size(data,3);
IMG(:,:,n)=fftshift(fft2(data(:,:,n)));
end




if no_rescale==1
  if nargin<3
     noise=abs(std(reshape(IMG(1:20,1:20,:),[],size(IMG,3))));
  end
 for n=1:size(IMG,3)
   IMG(:,:,n)=IMG(:,:,n)./(noise(n)./noise(1));
 end
end

RSOS=sqrt(sum(abs(IMG).^2,3));

% PHASE=sum(angle(IMG,3));




