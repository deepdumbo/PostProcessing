function [ output ] = remove_kspace_from_other_respiratory_stage(table_input_resp , j,  table_input_ecg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
            
            %% on recupere la table de ecg
            output=table_input_ecg;
            
%             q=find(output(:)>0);
            
            %% on recupere les indices qui ne correspondent pas au bon moment
            q=find(table_input_resp(:)~=j);
            
            %% on les supprime
            output(q)=0;
            
%             q=find(output(:)>0);

end

