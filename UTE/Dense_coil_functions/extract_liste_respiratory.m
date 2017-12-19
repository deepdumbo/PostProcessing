function [ table ] = extract_liste_respiratory(input , debut, fin, nb_phase_respiratoire )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% on exclut le debut et la fin de l'acquisition
navigator=input(debut:fin);

table.respiration.raw=single(zeros(size(navigator,1),nb_phase_respiratoire));
table.respiration.affichage=single(zeros(size(navigator,1),nb_phase_respiratoire));
table.respiration.liste=single(zeros(size(navigator,1),1));

% table.respiration_quartile.raw=zeros(size(navigator,1),nb_phase_respiratoire);
% table.respiration_quartile.affichage=zeros(size(navigator,1),nb_phase_respiratoire);
% table.respiration_quartile.liste=zeros(size(navigator,1),1);

for j=1:nb_phase_respiratoire
    
    de_ici=(j-1)/nb_phase_respiratoire;
    a_ici=(j)/nb_phase_respiratoire;
    
    for i=1:1:size(navigator,1)
        
        if ((navigator(i,1)>=de_ici) && (navigator(i,1)<a_ici))
            
            table.respiration.raw(i,j)=1;
            table.respiration.affichage(i,j)=1;
            table.respiration.liste(i,1)=j;
        else
            table.respiration.raw(i,j)=0;
            table.respiration.affichage(i,j)=-1;
            
        end
        
    end
    
%     de_ici2=quantile(navigator,de_ici);
%     a_ici2=quantile(navigator,a_ici);
%     
%     
%     for i=1:1:size(navigator,1)
%         
%         if ((navigator(i,1)>=de_ici2) && (navigator(i,1)<a_ici2))
%             
%             table.respiration_quartile.raw(i,j)=1;
%             table.respiration_quartile.affichage(i,j)=1;
%             table.respiration_quartile.liste(i,1)=j;
%         else
%             table.respiration_quartile.raw(i,j)=0;
%             table.respiration_quartile.affichage(i,j)=-1;
%         end
%         
%     end
    
    str_smg=sprintf('%d %2.3f %2.3f  %d %2.2f',j , de_ici, a_ici, ...
        sum(table.respiration.raw(:,j)), sum(table.respiration.raw(:,j))/size(table.respiration.raw,1) );  disp(str_smg);
%     str_smg=sprintf('%d %2.3f %2.3f  %d %2.2f',j , de_ici2, a_ici2, ...
%         sum(table.respiration_quartile.raw(:,j)), sum(table.respiration_quartile.raw(:,j))/size(table.respiration_quartile.raw,1) );  disp(str_smg);
    
end



table.respiration.liste_ready=single(zeros(size(input,1),1));
table.respiration.liste_ready(debut:fin)=table.respiration.liste;

% table.respiration_quartile.liste_ready=single(zeros(size(input,1),1));
% table.respiration_quartile.liste_ready(debut:fin)=table.respiration_quartile.liste;




end

