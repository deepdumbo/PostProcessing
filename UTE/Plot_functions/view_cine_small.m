function [  ] = view_cine( input1, input2,  number_of_output_cardiac_phases , numero1, numero2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


close(figure(8))
figure(8)

taille=3;

indice=70;

while(1)
    
    for ii=1:1:number_of_output_cardiac_phases
                
                      
            kk=35;
            subplot(2,taille, kk-34 +(1-1)*taille);
            temp=squeeze(input1(kk,:,:,ii,numero1));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=45;
            subplot(2,taille, kk-44+1 +(1-1)*taille);
            temp=squeeze(input1(:,kk,:,ii,numero1));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=55;
            subplot(2,taille, kk-54+2 +(1-1)*taille);
            temp=squeeze(input1(:,kk,:,ii,numero1));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=35;
            subplot(2,taille, kk-34 +(2-1)*taille);
            temp=squeeze(input2(kk,:,:,ii,numero2));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=45;
            subplot(2,taille, kk-44+1 +(2-1)*taille);
            temp=squeeze(input2(:,kk,:,ii,numero2));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            
            kk=55;
            subplot(2,taille, kk-54+2 +(2-1)*taille);
            temp=squeeze(input2(:,kk,:,ii,numero2));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
            %             end
            
       
        pause(0.01)
    end
end


end

