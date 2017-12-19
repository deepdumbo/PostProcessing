folder='/home/valeryo/mount/Imagerie/DICOM_DATA/Bruker_7617/12/';

ex = read_bru_experiment(folder);

%matrix
nX=ex.method.PVM_EncMatrix(1);
nY=ex.method.PVM_EncMatrix(2);
nZ=ex.method.PVM_EncMatrix(3);%condition si 3D

%bytorder
order=ex.acqp.BYTORDA;

%frequency
freq=ex.method.PVM_FrqRef(1);

%coils
nCoils=ex.method.PVM_EncNReceivers;
%encoding steps
ex.method.PVM_EncSteps1;
ex.method.PVM_EncSteps2;

nS=ex.method.PVM_SPackArrNSlices;%slices
nPS=ex.method.PVM_NSPacks;%pack of slices

nR=ex.method.PVM_NRepetitions;%repetitions

nE=ex.method.PVM_NEchoImages;%nombre d'echos
nEC=ex.method.PVM_EvolutionCycles;%nombre d'evolution ou temps d'inversions

%rare factor
nRF=ex.acqp.ACQ_rare_factor;

%aliasing
ex.method.PVM_AntiAlias;

%partial fourier encoding
ex.method.PVM_EncPft(1); %en read
ex.method.PVM_EncPft(2); %en phase

%object order scheme
ex.method.PVM_ObjOrderScheme;

%spatial phase
ex.acqp.ACQ_spatial_phase_1;
ex.acqp.ACQ_spatial_phase_2; %si 3D

reshape(ex.fid, [nX,nCoils,nY,nZ]);