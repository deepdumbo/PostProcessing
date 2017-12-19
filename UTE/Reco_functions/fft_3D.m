function [ out ] = fft_3D( in )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Reconstruct in x 
K = fftshift(fft(fftshift(in,1),[],1),1);

% Reconstruct in y 
temp = fftshift(fft(fftshift(K,2),[],2),2);

% Reconstruct in z
out = fftshift(fft(fftshift(temp,3),[],3),3);

end

