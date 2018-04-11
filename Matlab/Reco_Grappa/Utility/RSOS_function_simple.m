function [RSOS]=RSOS_function_simple(data)
     if nargin<2; no_rescale=0; end
%
% Convert k-space to RSOS
%

for n=1:size(data,3);
IMG(:,:,n)=fftshift(fft2(data(:,:,n)));
end


RSOS=sqrt(sum(abs(IMG).^2,3));

