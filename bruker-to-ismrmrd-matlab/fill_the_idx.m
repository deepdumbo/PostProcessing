function [ idx, flag ] = fill_the_idx( header , ex)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% get some info with readable name

% [ nX, nY, nZ ] = get_dimensions_acq( ex );
% nous voulons recuperer les dimensionstotales de la matrice d'encodage
[ nX, nY, nZ ] = get_dimensions( ex );

[ readout, E1, E2 ] = get_encoding_size( ex , nZ );

number_of_repetitions=ex.method.PVM_NRepetitions;
number_of_slices=ex.method.PVM_SPackArrNSlices;
number_of_averages=ex.method.PVM_NAverages;
number_of_stacks=ex.method.PVM_NSPacks;

if(isfield(ex.method,'PVM_NEchoImages'))
    number_of_echos=ex.method.PVM_NEchoImages;
else
    number_of_echos=1;
end

if(isfield(ex.method,'PVM_DwNDiffExp'))
    number_of_diffusion=ex.method.PVM_DwNDiffExp;
else
    number_of_diffusion=1;
end



clear idx

count=0;

for rep = 1:number_of_repetitions
    for  e2 = 1:1:E2       
          for  e1 = 1:1:E1
              for  ne = 1:1:number_of_echos
                for nd=1:1:number_of_diffusion
                    for  ns = 1:1:number_of_slices
                        
                        count = count +1;
                        
                        idx.kspace_encode_step_1(count)=ex.method.PVM_EncSteps1(e1) + round(nY/2);                       
                        
                        if (E2>1)
                            idx.kspace_encode_step_2(count)=ex.method.PVM_EncSteps2(e2)+ round(nZ/2);
                        else
                            idx.kspace_encode_step_2(count)=1-1;
                        end
                        
                        %                 if(
                        %                 idx.parallel_and_calib(count)=
                        
                        idx.slice(count)=ns-1;
                        idx.repetition(count)=rep-1;                       
                        idx.set(count)=nd-1;                                           
                        idx.contrast(count)=ne-1;
                       
                        
                        %                     if (e1==1)
%                        str_msg=sprintf('count %d slice %d e1 %d e2 %d ne %d  rep %d  ', count, idx.slice(count), idx.kspace_encode_step_1(count), idx.kspace_encode_step_2(count), idx.contrast(count), idx.repetition(count)); disp( str_msg);
                        %                     end
                        
                        if (e1==1)
                            flag.first_in_encoding_step1(count)=1;
                        else
                            flag.first_in_encoding_step1(count)=0;
                        end
                        
                        if (e1==E1)
                            flag.last_in_encoding_step1(count)=1;
                        else
                            flag.last_in_encoding_step1(count)=0;
                        end
                        
                        
                        if (e1==1)
                            flag.first_in_contrast(count)=1;
                        else
                            flag.first_in_contrast(count)=0;
                        end
                        
                        if (e1==E1)
                            flag.last_in_contrast(count)=1;
                        else
                            flag.last_in_contrast(count)=0;
                        end
                        
                        %                 if (e1==1)
                        %                 flag.first_in_encoding(count)=1;
                        %                 else
                        %                 flag.first_in_encoding(count)=0;
                        %                 end
                    end
                end
            end
        end
    end
end



end

