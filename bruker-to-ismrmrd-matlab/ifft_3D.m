function [ out ] = ifft_3D( in )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Reconstruct in x 
K = fftshift(ifft(fftshift(in,1),[],1),1);

% Reconstruct in y 
temp=fftshift(ifft(fftshift(K,2),[],2),2);

% Reconstruct in z 
out=fftshift(ifft(fftshift(temp,3),[],3),3);

% je reduits
% [R0 E1 E2 CHA N S SLC ] - > [ R0 E1 E2 CHA]

% remplacer N ou S par ECHO ou contrast


end

