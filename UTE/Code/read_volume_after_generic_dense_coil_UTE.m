clear all 
close all

A=load('/home/valery/TempoFigures/image_resp_1_cardiac_20_filter_N.mat');

volume_1=A.volume; clear A 

size(volume_1)

A=load('/home/valery/TempoFigures/image_resp_3_cardiac_20_filter_N.mat');

volume_2=A.volume; clear A 

size(volume_2)

A=load('/home/valery/TempoFigures/image_resp_3_cardiac_20_filter_Y.mat');

volume_3=A.volume; clear A 

size(volume_3)












number_of_output_cardiac_phases=20;
number_of_output_respiratory_phases=3;

view_cine( volume_2, volume_3, number_of_output_cardiac_phases , 1, 1);



ii=2
indice=30;

while(1)
    for jj=1:1:number_of_output_respiratory_phases
        for kk=indice:10:130
            
            subplot(number_of_output_respiratory_phases,4, kk/10-2);
            temp=squeeze(volume_4(kk,:,:,ii,jj));
            imagesc(imrotate(temp,90));  title( int2str(kk)); colormap(gray);
        end
        
        pause(0.4)
    end
    
end


view_from_base_to_apex( volume_2(:,:,:,1) , volume_4(:,:,:,1) );
 
view_from_sagital( volume_2(:,:,:,1,1) , volume_4(:,:,:,1) );

view_cine( volume_2, volume_3, number_of_output_cardiac_phases , 1, 1);







 
 

%% partie moco 
close(figure(11))
figure(11)

ii=1;

 for jj=1:1:number_of_output_respiratory_phases
     
     subplot(2,3,jj);
     temp=squeeze(volume_2(110,:,:,ii,jj));
     imagesc(temp);  title( int2str(kk)); colormap(gray);
     
     subplot(2,3,jj+3);
     imagesc(temp(85:95,:));  title( int2str(kk)); colormap(gray);  
     
     lala(:,jj)=mean(temp(80:100,:),1);
      
 end

 

x=[62 , 58 , 53 ];

deplacement=x-x(1);


volume_2_new=zeros(size(volume_2));

 for jj=1:1:number_of_output_respiratory_phases

volume_2_new(:,:,11+deplacement(jj):end-10+deplacement(jj),ii,jj)=volume_2(:,:,11:end-10,ii,jj);

 end
 
 
 
%% partie moco 
close(figure(11))
figure(11)

ii=1;

 for jj=1:1:number_of_output_respiratory_phases
     
     subplot(2,3,jj);
     temp=squeeze(volume_2_new(110,:,:,ii,jj));
     imagesc(temp);  title( int2str(kk)); colormap(gray);
     
     subplot(2,3,jj+3);
     imagesc(temp(85:95,:));  title( int2str(kk)); colormap(gray);  
     
     lala_corrige(:,jj)=mean(temp(80:100,:),1);
      
 end
 
 
 
close(figure(12))
figure(12)

 subplot(1,2,1);
 
for jj=1:1:number_of_output_respiratory_phases
 
 plot(lala(:,jj)); hold on;
end
 
 subplot(1,2,2);
for jj=1:1:number_of_output_respiratory_phases
 
 plot(lala_corrige(:,jj)); hold on;
end

volume_mean=mean(volume_2_new(:,:,:,ii,:),5);


view_from_base_to_apex( volume_2(:,:,:,ii,1) , volume_mean(:,:,:,ii) );





A=load('/home/valery/TempoFigures/image_resp_3_cardiac_1_filter_N.mat');

volume_resp_3=A.volume; clear A 

size(volume_resp_3)



A=load('/home/valery/TempoFigures/image_resp_5_cardiac_1_filter_N.mat');

volume_resp_5=A.volume; clear A 

size(volume_resp_5)


A=load('/home/valery/TempoFigures/image_resp_7_cardiac_1_filter_N.mat');

volume_resp_7=A.volume; clear A 

size(volume_resp_7)



close(figure(13))
figure(13)

for jj=1:1:size(volume_resp_3,5)
    
     temp=squeeze(volume_resp_3(110,:,:,1,jj));     
     subplot(1,size(volume_resp_3,5),jj); imagesc(temp); colormap(gray);
    
end

close(figure(14))
figure(14)

for jj=1:1:size(volume_resp_5,5)
    
     temp=squeeze(volume_resp_5(110,:,:,1,jj));     
     subplot(1,size(volume_resp_5,5),jj); imagesc(temp); colormap(gray);
    
end


close(figure(15))
figure(15)

for jj=1:1:size(volume_resp_7,5)
    
     temp=squeeze(volume_resp_7(110,21:148,1:128,1,jj));     
     subplot(1,size(volume_resp_7,5),jj); imagesc(temp); colormap(gray);
    
end
