function [ idx, flag ] = fill_the_idx_mp2rage( header , ex)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




% sy=ex.method.PVM_Matrix(2);
% sz=ex.method.PVM_Matrix(3);
% GradPhaseVector=ex.method.GradPhaseVector;
% GradSliceVector=ex.method.GradSliceVector;
% 
% GradPhaseVector=GradPhaseVector*sy/2;
% GradPhaseVector=round(GradPhaseVector-min(GradPhaseVector)+1);
% 
% GradSliceVector=GradSliceVector*sz/2;
% GradSliceVector=round(GradSliceVector-min(GradSliceVector)+1);
% % mask
% mask=zeros(sy,sz);
% 
% mask(sub2ind(size(mask),GradPhaseVector(:),GradSliceVector(:)))=1;
% figure;imshow(mask);

%% get some info with readable name

[ nX, nY, nZ ] = get_dimensions_acq( ex );

[ readout, E1, E2 ] = get_encoding_size( ex , nZ );

number_of_repetitions=ex.method.PVM_NRepetitions;
number_of_slices=ex.method.PVM_SPackArrNSlices;
number_of_averages=ex.method.PVM_NAverages;
number_of_echos=ex.method.PVM_NEchoImages;
number_of_stacks=ex.method.PVM_NSPacks;



GradPhaseVector=ex.method.GradPhaseVector*E1/2;
GradPhaseVector=round(GradPhaseVector-min(GradPhaseVector)+1);

GradSliceVector=ex.method.GradSliceVector*E2/2;
GradSliceVector=round(GradSliceVector-min(GradSliceVector)+1);

GradPhaseVector2=reshape(GradPhaseVector,ex.acqp.ACQ_rare_factor,[]);
GradPhaseVector2=reshape([GradPhaseVector2 GradPhaseVector2],[],1);

GradSliceVector2=reshape(GradSliceVector,ex.acqp.ACQ_rare_factor,[]);

if (strcmp(ex.method.EnableCS,'Fully') ||  strcmp(ex.method.EnableCS,'No'))
    GradSliceVector2=reshape([GradSliceVector2; GradSliceVector2],[],1);
else
    GradSliceVector2=reshape([GradSliceVector2 GradSliceVector2],[],1);
end

temp1=ones(ex.acqp.ACQ_size(2),1);
temp2=2*ones(ex.acqp.ACQ_size(2),1);

temp1=reshape(temp1,ex.acqp.ACQ_rare_factor,[]);
temp2=reshape(temp2,ex.acqp.ACQ_rare_factor,[]);

ContrastVector=reshape([temp1;temp2],[],1);

% mask
mask=zeros(E1,E2);

mask(sub2ind(size(mask),GradPhaseVector(:),GradSliceVector(:)))=1;
figure;imshow(mask);

% count=0;
% for  ll = 1:1:size(GradPhaseVector2,1)
% 
%  count=count+1;
%  
%  idx.kspace_encode_step_1(count)=GradPhaseVector2(ll)-1;
%  idx.kspace_encode_step_2(count)=GradSliceVector2(ll)-1;
% 
%  idx.slice(count)=0;
%  idx.repetition(count)=0;                    
%  idx.contrast(count)=ContrastVector(ll)-1;
%  
%  str_msg=sprintf('count %d e1 %d e2 %d  echo %d   ', count,  idx.kspace_encode_step_1(count), idx.kspace_encode_step_2(count), idx.contrast(count)); disp( str_msg);
%   
% end;


%% test aurel

disp('fin gestion encodage');

count = 0;
counter2=0;
for k=1:size(GradPhaseVector,1)/ex.acqp.ACQ_rare_factor
        for i=1:ex.acqp.ACQ_rare_factor
             count=count+1;
             counter2=counter2+1;
             
             idx.kspace_encode_step_1(count)=GradPhaseVector(counter2)-1;
             idx.kspace_encode_step_2(count)=GradSliceVector(counter2)-1;
             idx.slice(count)=0;
             idx.repetition(count)=0;
             idx.contrast(count)=0;
             
             idx.kspace_encode_step_1(count+ex.acqp.ACQ_rare_factor)=GradPhaseVector(counter2)-1;
             idx.kspace_encode_step_2(count+ex.acqp.ACQ_rare_factor)=GradSliceVector(counter2)-1;
             idx.slice(count+ex.acqp.ACQ_rare_factor)=0;
             idx.repetition(count+ex.acqp.ACQ_rare_factor)=0;
             idx.contrast(count+ex.acqp.ACQ_rare_factor)=1;
        end
count=count+ex.acqp.ACQ_rare_factor;
end

disp('fin gestion encodage');


end

