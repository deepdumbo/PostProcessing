function [  ] = view_cine( input1, input2,  number_of_output_cardiac_phases , numero1, numero2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


close(figure(8))
figure(8)

taille=3;

indice=70;

while(1)
    
    for ii=1:1:number_of_output_cardiac_phases
                
                      
            kk=65;
            subplot(2,taille, kk-64 +(1-1)*taille);
            temp=squeeze(input1(kk,:,:,ii,numero1));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=75;
            subplot(2,taille, kk-74+1 +(1-1)*taille);
            temp=squeeze(input1(:,kk,:,ii,numero1));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=90;
            subplot(2,taille, kk-89+2 +(1-1)*taille);
            temp=squeeze(input1(:,kk,:,ii,numero1));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=65;
            subplot(2,taille, kk-64 +(2-1)*taille);
            temp=squeeze(input2(kk,:,:,ii,numero2));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=75;
            subplot(2,taille, kk-74+1 +(2-1)*taille);
            temp=squeeze(input2(:,kk,:,ii,numero2));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=90;
            subplot(2,taille, kk-89+2 +(2-1)*taille);
            temp=squeeze(input2(:,kk,:,ii,numero2));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            %             end
            
       
        pause(0.03)
    end
end


end

