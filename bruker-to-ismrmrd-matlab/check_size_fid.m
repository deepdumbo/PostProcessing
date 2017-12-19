function [  ] = check_size_fid( ex  , data_for_acqp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid_size=prod(size(data_for_acqp));

% acqp_size=prod(ex.acqp.ACQ_size)/2*ex.method.PVM_EncNReceivers;
% 
% if (fid_size~=acqp_size)
%    
%     fprintf('%d %d \n', fid_size, acqp_size);
%     disp('PROBLEME : la taille de la matrice fid ne correspond pas à ACQ_size');
% end  
%     
% acqp_size=prod(ex.acqp.ACQ_size)/2*ex.method.PVM_EncNReceivers*ex.method.PVM_SPackArrNSlices;
%     
% if (fid_size~=acqp_size)     
%     
%     fprintf('%d %d \n', fid_size, acqp_size);
%     disp('PROBLEME : la taille de la matrice fid ne correspond pas à ACQ_size');
%     
% end
% 
% acqp_size=prod(ex.acqp.ACQ_size)/2*ex.method.PVM_EncNReceivers*ex.method.PVM_SPackArrNSlices*ex.method.PVM_NRepetitions;
%     
% if (fid_size~=acqp_size)     
%     
%     fprintf('%d %d \n', fid_size, acqp_size);
%     disp('PROBLEME : la taille de la matrice fid ne correspond pas à ACQ_size');
%     
% end



[ number_of_channels , number_of_repetitions , number_of_slices , number_of_echos , number_of_diffusion] = find_number_of( ex );

acqp_size=prod(ex.acqp.ACQ_size)/2*number_of_channels*number_of_slices*number_of_repetitions*number_of_diffusion;
    
if (fid_size~=acqp_size)     
    
    fprintf('%d %d \n', fid_size, acqp_size);
    disp('PROBLEME : la taille de la matrice fid ne correspond pas à ACQ_size');
    
else
    disp('La taille de la matrice fid correspond à ACQ_size');
    
end




end

