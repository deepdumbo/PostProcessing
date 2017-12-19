function [ output_args ] = view_from_base_to_apex( input1 , input2 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



figure()

for ii=50:2:size(input1,2)-50
    
       
        subplot(1,2,1)
        temp=squeeze(input1(:,ii,:));
        imagesc(imrotate(temp,90));  title( int2str(ii)); colormap(gray);
        subplot(1,2,2)
        temp=squeeze(input2(:,ii,:));
        imagesc(imrotate(temp,90));  title( int2str(ii)); colormap(gray);
    
    pause(0.5)
end

close all
end

