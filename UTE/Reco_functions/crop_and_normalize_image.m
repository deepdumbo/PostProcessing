function [ output ] = crop_and_normalize_image( input, osfcrop , Ncrop , factor )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


            %% crop image  
            image_tempo=input(round((round(osfcrop*Ncrop)-Ncrop)/2+1):round(Ncrop+(round(osfcrop*Ncrop)-Ncrop)/2),round((round(osfcrop*Ncrop)-Ncrop)/2+1):round(Ncrop+(round(osfcrop*Ncrop)-Ncrop)/2),round((round(osfcrop*Ncrop)-Ncrop)/2+1):round(Ncrop+(round(osfcrop*Ncrop)-Ncrop)/2));
            
            %% IMAGE OUPUT
            output(:,:,:)=(abs(image_tempo)-min(min(min(abs(image_tempo)))))* factor/(max(max(max(abs(image_tempo))))-min(min(min(abs(image_tempo)))));


end

