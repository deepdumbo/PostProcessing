function [ output ] = start_list_ecg_in_systole( input , debut_change, nb_phase_ecg)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

output=zeros(size(input));


for jj=debut_change:nb_phase_ecg
    q=find(input(:,1)==jj);
    output(q)=jj-(debut_change-1);
%     jj-(debut_change-1)
end

for jj=1:debut_change-1
    q=find(input(:,1)==jj);
    output(q)=jj+nb_phase_ecg-(debut_change-1);
%     jj+nb_phase_ecg-(debut_change-1)
end

end

