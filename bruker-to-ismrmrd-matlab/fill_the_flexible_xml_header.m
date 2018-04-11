function [ header ] = fill_the_flexible_xml_header( ex )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [];

[ nX, nY, nZ ] = get_dimensions( ex );

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = ex.method.PVM_FrqRef(1)*1e6; % 9.4T


header.measurementInformation.measurementID='2';
header.measurementInformation.patientPosition='HFS';
header.measurementInformation.protocolName='Flash';
header.measurementInformation.frameOfReferenceUID='lalala'
header.measurementInformation.measurementDependency.dependencyType='Noise';
header.measurementInformation.measurementDependency.measurementID='1';



% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = ex.acqp.ORIGIN{1};
header.acquisitionSystemInformation.systemModel = ex.acqp.ACQ_station{1};
header.acquisitionSystemInformation.systemFieldStrength_T= ex.acqp.BF1/42.577;
header.acquisitionSystemInformation.receiverChannels = ex.method.PVM_EncNReceivers;
header.acquisitionSystemInformation.relativeReceiverNoiseBandwidth=1;
header.acquisitionSystemInformation.institutionName =ex.acqp.ACQ_institution{1};
header.acquisitionSystemInformation.stationName =ex.acqp.ACQ_station{1};
 

fprintf('-------------------------------------------------------------- \n');
fprintf('%s\n', ex.method.Method{1,1});
fprintf('-------------------------------------------------------------- \n');
fprintf('%s \n',header.acquisitionSystemInformation.systemVendor);
fprintf('%s \n',header.acquisitionSystemInformation.systemModel);
fprintf('%s \n',header.acquisitionSystemInformation.systemFieldStrength_T);
fprintf('%s \n',header.acquisitionSystemInformation.institutionName);
fprintf('%s \n',header.acquisitionSystemInformation.stationName);
fprintf('-------------------------------------------------------------- \n');


% The Encoding (Required)
header.encoding.trajectory = 'cartesian';
header.encoding.encodedSpace.fieldOfView_mm.x = ex.acqp.ACQ_fov(1)*10;
header.encoding.encodedSpace.fieldOfView_mm.y = ex.acqp.ACQ_fov(2)*10;
header.encoding.encodedSpace.matrixSize.x = nX;
header.encoding.encodedSpace.matrixSize.y = nY;

if(nZ>1)
    header.encoding.encodedSpace.fieldOfView_mm.z = ex.acqp.ACQ_fov(3)*10;
    header.encoding.encodedSpace.matrixSize.z = nZ;
else
    header.encoding.encodedSpace.fieldOfView_mm.z = 1;
    header.encoding.encodedSpace.matrixSize.z = 1;
end


% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;

% header.encoding.reconSpace.matrixSize.x=header.encoding.encodedSpace.matrixSize.x*2;
% header.encoding.reconSpace.matrixSize.y=header.encoding.encodedSpace.matrixSize.y*2;
% if (header.encoding.encodedSpace.matrixSize.z~=1)
%     header.encoding.reconSpace.matrixSize.z=header.encoding.encodedSpace.matrixSize.z;
% end


if (ex.method.PVM_EncPft(1)~=1)
    disp('attention partial fourier suivant x')
end

if (ex.method.PVM_EncPft(2)~=1)
    disp('attention partial fourier suivant y')
    header.encoding.encodingLimits.kspace_encoding_step_1.minimum = ex.method.PVM_EncSteps1(1) + round(nY/2);
    header.encoding.encodingLimits.kspace_encoding_step_1.maximum =   ex.method.PVM_EncSteps1(end) + round(nY/2);
    header.encoding.encodingLimits.kspace_encoding_step_1.center =  round(nY/2)+1;
else
    header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_1.maximum = nY-1;
    header.encoding.encodingLimits.kspace_encoding_step_1.center = floor(nY/2);
end


% Encoding Limits
% header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0;
% header.encoding.encodingLimits.kspace_encoding_step_0.maximum = nX-1;
% header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(nX/2);


if(nZ>1)
    header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = nZ-1;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = floor(nZ/2);
    
else
    header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum =0;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = 0;
    
end

[ number_of_slices ] = get_number_of_slices( ex );

% if(number_of_slices>1)
%     header.encoding.encodingLimits.slice.minimum = 0;
%     header.encoding.encodingLimits.slice.maximum = nZ-1;
%     header.encoding.encodingLimits.slice.center = floor(nZ/2)-1;    
% end    
    
 if (number_of_slices>1)
    header.encoding.encodingLimits.slice.minimum = 0;
    header.encoding.encodingLimits.slice.maximum = ex.method.PVM_SPackArrNSlices-1;
    header.encoding.encodingLimits.slice.center = 0;
    fprintf('number of slices is %d \n', ex.method.PVM_SPackArrNSlices);
end

if(isfield(ex.method,'PVM_NRepetitions'))
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = ex.method.PVM_NRepetitions-1;
header.encoding.encodingLimits.repetition.center = 0;
fprintf('number of repetitions %d \n', ex.method.PVM_NRepetitions);
else 
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = 0;
header.encoding.encodingLimits.repetition.center = 0;    
end

% header.encoding.encodingLimits.average.minimum = 0;
% header.encoding.encodingLimits.average.maximum = ex.method.PVM_NAverages-1;
% header.encoding.encodingLimits.average.center = 0;


if(isfield(ex.method,'PVM_NEchoImages'))
header.encoding.encodingLimits.contrast.minimum = 0;
header.encoding.encodingLimits.contrast.maximum = ex.method.PVM_NEchoImages-1;
header.encoding.encodingLimits.contrast.center = 0;
fprintf('number of echos %d \n', ex.method.PVM_NEchoImages);
else    
header.encoding.encodingLimits.contrast.minimum = 0;
header.encoding.encodingLimits.contrast.maximum = 0;
header.encoding.encodingLimits.contrast.center = 0;
end

if(isfield(ex.method,'PVM_DwNDiffExp'))
header.encoding.encodingLimits.set.minimum = 0;
header.encoding.encodingLimits.set.maximum = ex.method.PVM_DwNDiffExp-1;
header.encoding.encodingLimits.set.center = 0;
fprintf('number of echos diffusion%d \n', ex.method.PVM_DwNDiffExp);
else    
header.encoding.encodingLimits.set.minimum = 0;
header.encoding.encodingLimits.set.maximum = 0;
header.encoding.encodingLimits.set.center = 0;
end


fprintf('-------------------------------------------------------------- \n');
fprintf('Number of encoding spaces : %d \n',nX);
fprintf('Encoding matrix size : %d %d %d\n', header.encoding.encodedSpace.matrixSize.x,  header.encoding.encodedSpace.matrixSize.y ,  header.encoding.encodedSpace.matrixSize.z);
fprintf('Encoding field_of_view : %d %d %d\n', header.encoding.encodedSpace.fieldOfView_mm.x, header.encoding.encodedSpace.fieldOfView_mm.y , header.encoding.encodedSpace.fieldOfView_mm.z);
fprintf('Recon matrix size : %d %d %d\n', header.encoding.reconSpace.matrixSize.x,  header.encoding.reconSpace.matrixSize.y ,  header.encoding.reconSpace.matrixSize.z);
fprintf('Recon field_of_view : %d %d %d\n', header.encoding.reconSpace.fieldOfView_mm.x, header.encoding.reconSpace.fieldOfView_mm.y , header.encoding.reconSpace.fieldOfView_mm.z);
fprintf('-------------------------------------------------------------- \n');

header.sequenceParameters.TR= ex.method.PVM_RepetitionTime;
header.sequenceParameters.TE= ex.acqp.ACQ_echo_time;

% header.sequenceParameters.TE(2)= 4;
% header.sequenceParameters.TE(3)= ex.method.PVM_EchoTime+2;
% header.sequenceParameters.TE(4)= ex.method.PVM_EchoTime+3;
header.sequenceParameters.flipAngle_deg=ex.acqp.ACQ_flip_angle;




end

