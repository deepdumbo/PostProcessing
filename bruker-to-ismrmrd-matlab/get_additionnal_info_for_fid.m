function [ readout, number_of_channels ] = get_additionnal_info_for_fid( ex )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%
readout=ex.method.PVM_EncMatrix(1);

%
number_of_channels=ex.method.PVM_EncNReceivers;

nombre_de_phases_suivant_acq_y=ex.acqp.ACQ_spatial_size_1;

nombre_de_phases_suivant_y=size(ex.method.PVM_EncSteps1,1);
nombre_de_phases_suivant_z=size(ex.method.PVM_EncSteps2,1);

% nombre_de_points_par_antennes=size(ex.fid,1)/number_of_channels;

% number_of_lines=size(ex.fid,1)/number_of_channels/readout;


% size(ex.fid,1)/nCoils/readout/E2;
% size(ex.fid,1)/nCoils/readout/E1;

end

