function [ output ] = remove_bad_peaks( data_peak_dectection , time_difference , repolarisation_time, Tr_, value )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


q=find(time_difference(:) < repolarisation_time);

compteur=0;

N=size(data_peak_dectection);



for i=N-1:-1:2;
    
    lala=find(q(:)==i);
    
    if (size(lala,1)>0)
                    
            str_smg= sprintf('%d %f %f %f \n', i, (data_peak_dectection(i,1)-data_peak_dectection(i-1,1))*Tr_ ,...
                (data_peak_dectection(i+1,1)-data_peak_dectection(i,1))*Tr_ ,  (data_peak_dectection(i+1,1)-data_peak_dectection(i-1,1))*Tr_ ); disp(str_smg);
            
            % double detection qui emepeche les de creer de grands intervalle
            if ((data_peak_dectection(i+1,1)-data_peak_dectection(i-1,1))*Tr_> value)
                
                % on ne jette pas les donnes
                compteur=compteur+1;
                nouveau_data_peak(compteur,1)=data_peak_dectection(i,1);
                nouveau_data_peak(compteur,2)=data_peak_dectection(i,2);
                
                 str_smg= sprintf('%d on ne sait pas on garde \n', i); disp(str_smg);
                 
            else
                
                % on jette les donnes
                 str_smg= sprintf('%d on jette \n', i); disp(str_smg);
            end     
        
    else
        compteur=compteur+1;
        nouveau_data_peak(compteur,1)=data_peak_dectection(i,1);
        nouveau_data_peak(compteur,2)=data_peak_dectection(i,2);
    end
    
end


output=nouveau_data_peak(end:-1:1,:);


end

