clear all


[status,id]= system('whoami');

str_user= id(1:end-1);

check_if_iam_using_the_ihuserver(str_user);

[ str_network_imagerie, str_network_perso ] = get_network_name( str_user );

% filename='/home/valery/generic_spirit_random_cine.h5'
filename=['/home/',str_user,str_network_imagerie , '/DICOM_DATA/2017-10-31_SMS_Brain_Aurelien/FIDSMS/meas_MID00110_FID51718_cmrr_1_7iso_GP2_MB2_36s_50sur.h5'];

% filename='/home/valery/testdata.h5';
% filename='/home/valery/Reseau/Imagerie/For_Valery/testdata.h5'

if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them


hdr = ismrmrd.xml.deserialize(dset.readxml);

%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

% Number of slices, coils, repetitions, contrasts etc.
% We have to wrap the following in a try/catch because a valid xml header may
% not have an entry for some of the parameters

try
    number_of_slices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    number_of_slices = 1;
end

try
    number_of_coils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    number_of_coils = 1;
end

try
    number_of_repetitions = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    number_of_repetitions = 1;
end

try
    number_of_contrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 + 1;
catch
    number_of_contrasts = 1;
end

try
    number_of_sets = hdr.encoding.encodingLimits.set.maximum + 1 + 1;
catch
    number_of_sets = 1;
end

D = dset.readAcquisition();

% Note: can select a single acquisition or header from the block, e.g.
%  acq = D.select(1);
%  hdr = D.head.select(1);
% or you can work with them all together





%% Ignore noise scans
% TODO add a pre-whitening example
% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
% isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
% firstScan = find(isNoise==0,1,'first');
% if firstScan > 1
%     noise = D.select(1:firstScan-1);
% else
%     noise = [];
% end

maximum_scan_number=3637;
meas  = D.select(1:maximum_scan_number);
% clear D;


% ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1               =  1,
% ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1                =  2,
% ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2               =  3,
% ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2                =  4,
% ISMRMRD_ACQ_FIRST_IN_AVERAGE                    =  5,
% ISMRMRD_ACQ_LAST_IN_AVERAGE                     =  6,
% ISMRMRD_ACQ_FIRST_IN_SLICE                      =  7,
% ISMRMRD_ACQ_LAST_IN_SLICE                       =  8,
% ISMRMRD_ACQ_FIRST_IN_CONTRAST                   =  9,
% ISMRMRD_ACQ_LAST_IN_CONTRAST                    = 10,
% ISMRMRD_ACQ_FIRST_IN_PHASE                      = 11,
% ISMRMRD_ACQ_LAST_IN_PHASE                       = 12,
% ISMRMRD_ACQ_FIRST_IN_REPETITION                 = 13,
% ISMRMRD_ACQ_LAST_IN_REPETITION                  = 14,
% ISMRMRD_ACQ_FIRST_IN_SET                        = 15,
% ISMRMRD_ACQ_LAST_IN_SET                         = 16,
% ISMRMRD_ACQ_FIRST_IN_SEGMENT                    = 17,
% ISMRMRD_ACQ_LAST_IN_SEGMENT                     = 18,
% ISMRMRD_ACQ_IS_NOISE_MEASUREMENT                = 19,
% ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION             = 20,
% ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING = 21,
% ISMRMRD_ACQ_IS_REVERSE                          = 22,
% ISMRMRD_ACQ_IS_NAVIGATION_DATA                  = 23,
% ISMRMRD_ACQ_IS_PHASECORR_DATA                   = 24,
% ISMRMRD_ACQ_LAST_IN_MEASUREMENT                 = 25,
% ISMRMRD_ACQ_IS_HPFEEDBACK_DATA                  = 26,
% ISMRMRD_ACQ_IS_DUMMYSCAN_DATA                   = 27,
% ISMRMRD_ACQ_IS_RTFEEDBACK_DATA                  = 28,
% ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA   = 29,
% 
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_ENCODE_STEP1')   );
% str_msg=sprintf(' ACQ_FIRST_IN_ENCODE_STEP1 , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_ENCODE_STEP1')   );
% str_msg=sprintf(' ACQ_LAST_IN_ENCODE_STEP1 , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_ENCODE_STEP2')   );
% str_msg=sprintf(' ACQ_FIRST_IN_ENCODE_STEP2 , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_ENCODE_STEP2')   );
% str_msg=sprintf(' ACQ_LAST_IN_ENCODE_STEP2 , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_AVERAGE')   ); 
% str_msg=sprintf(' ACQ_FIRST_IN_AVERAGE , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_AVERAGE')   );
% str_msg=sprintf(' ACQ_LAST_IN_AVERAGE , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_SLICE')   );
% str_msg=sprintf(' ACQ_FIRST_IN_SLICE , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_SLICE')   );
% str_msg=sprintf(' ACQ_LAST_IN_SLICE , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_CONTRAST')   );
% str_msg=sprintf(' ACQ_FIRST_IN_CONTRAST , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_CONTRAST')   );
% str_msg=sprintf(' ACQ_LAST_IN_CONTRAST , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_PHASE')   );
% str_msg=sprintf(' ACQ_FIRST_IN_PHASE , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_PHASE')   );
% str_msg=sprintf(' ACQ_LAST_IN_PHASE , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_REPETITION')   );
% str_msg=sprintf(' ACQ_FIRST_IN_REPETITION , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_REPETITION')   );
% str_msg=sprintf(' ACQ_LAST_IN_REPETITION , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_SET')   );
% str_msg=sprintf(' ACQ_FIRST_IN_SET , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_SET')   );
% str_msg=sprintf(' ACQ_LAST_IN_SET , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_FIRST_IN_SEGMENT')   );
% str_msg=sprintf(' ACQ_FIRST_IN_SEGMENT , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_LAST_IN_SEGMENT')   );
% str_msg=sprintf(' ACQ_LAST_IN_SEGMENT , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')   );
% str_msg=sprintf(' ACQ_IS_PARALLEL_CALIBRATION , %d ',  size(acq_test,2)); disp(str_msg);
% acq_test   = find(  meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING')   );
% str_msg=sprintf(' ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING , %d ',  size(acq_test,2)); disp(str_msg);




size(meas)


output_filename = '/home/valery/DICOM/testdata_no_MB.h5';

delete(output_filename)
% Create an empty ismrmrd dataset
if exist(output_filename,'file')
    error(['File ' output_filename ' already exists.  Please remove first'])
end

dset = ismrmrd.Dataset(output_filename);

dset.appendAcquisition(meas);

dset = ismrmrd.Dataset(output_filename);

xmlstring = ismrmrd.xml.serialize(hdr);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();


% 
% taille=size(meas.head.idx.kspace_encode_step_1,2)
% 
% vec=meas.head.idx.kspace_encode_step_1;
% 
% value_max=max(vec);
% 
% encoding_matrix=zeros(value_max+1,1);
% 
% 
% for i=1:1:taille
%     
%     if (meas.head.idx.set(i)==0)
%         
%         encoding_matrix(meas.head.idx.kspace_encode_step_1(i)+1,1)=encoding_matrix(meas.head.idx.kspace_encode_step_1(i)+1,1)+1;
%     end
% end
% 
% 
% plot(encoding_matrix)