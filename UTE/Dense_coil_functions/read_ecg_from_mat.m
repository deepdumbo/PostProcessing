function [ output ] = read_ecg_from_mat( filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ecg_read=load(filename);

ecg_tempo(:,:)=single(ecg_read.DonneesTableau(1:2:end,:));

[ output ] = cut_ecg_at_the_end_of_the_sequence( ecg_tempo );

output=single(output);

end

