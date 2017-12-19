

clear all
close all

number_of_output_respiratory_phases=3;
number_of_output_cardiac_phases=20;

% stop à 2 method , 2 resp , 2 cardiac

for method=1:1
    
    for j=1:1:number_of_output_respiratory_phases
        for i=1:number_of_output_cardiac_phases
            
            
            str_i=sprintf('%d',i );
            str_j=sprintf('%d',j );
            
            fprintf('method: %d resp: %d cardiac: %d\n',method, j, i );
            
            if (method==1)
                str_method='adapt';
            elseif (method==2)
                str_method='inati';
            else
                str_method='';
            end
            
            filename_img_cgsense=['/home/valery/Tempo/' , 'img_cgsense_', str_j , '_', str_i , '_',str_method , '.mat'];
            filename_img_comb=['/home/valery/Tempo/' , 'img_comb_', str_j , '_', str_i , '_',str_method ,  '.mat'];
            
            A=  load(filename_img_comb);
            B=  load(filename_img_cgsense);
            
            if (method==1)
                
                img_comb_adapt(:,:,:,i,j)=A.img_comb;
                img_cgsense_adapt (:,:,:,i,j)=B.img_cgsense;
                
%             elseif (method==2 && j==1)
%                 
%                 img_comb_inati(:,:,:,i,j)=A.img_comb;
%                 img_cgsense_inati(:,:,:,i,j)=B.img_cgsense;
                
            end            
        end
    end
end


addpath('/home/valery/Dev/cute/UTE/Code/');


view_cine( abs(img_comb_adapt), abs(img_cgsense_adapt), number_of_output_cardiac_phases , 1, 1);

view_cine( abs(img_cgsense_adapt), abs(img_cgsense_inati), number_of_output_cardiac_phases , 1, 1);

