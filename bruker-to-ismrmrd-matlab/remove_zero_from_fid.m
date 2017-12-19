function [ data_for_acqp ] = remove_zero_from_fid( ex )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[readout, number_of_channels]=get_additionnal_info_for_fid(ex);

ADC_size=128;

sizeR2=128*ceil(readout*number_of_channels/128);

acqp_begin=reshape(ex.fid, sizeR2,[]);

acqp_end=acqp_begin(1:(readout*number_of_channels),:);

number_of_lines=size(acqp_end,2);

fprintf('nombre_de_points_par_antennes %d \n', prod(size(acqp_end))/number_of_channels) ; 
fprintf('nombre_de_lignes %d \n', number_of_lines);

data_for_acqp=reshape(acqp_end,[readout,number_of_channels,number_of_lines]);

check_size_fid(ex, data_for_acqp);








end

