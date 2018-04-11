function [ idx, flag ] = fill_the_idx( header , ex)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% get some info with readable name

% nous voulons recuperer les dimensions de la matrice d'encodage

[ nX, nY, nZ ] = get_dimensions( ex );

[ readout, E1, E2 ] = get_encoding_size( ex , nZ );

[ number_of_repetitions ] = get_number_of_repetitions( ex );

[ number_of_slices ] = get_number_of_slices( ex );

[ number_of_averages ] = get_number_of_averages( ex );

[ number_of_stacks ] = get_number_of_stacks( ex );

[ number_of_echos ] = get_number_of_echos( ex );

[ number_of_diffusion ] = get_number_of_diffusion( ex );

clear idx

count=0;

for rep = 1:number_of_repetitions
    for  e2 = 1:1:E2
        for  e1 = 1:1:E1
            for  ne = 1:1:number_of_echos
                for nd=1:1:number_of_diffusion
                    for  ns = 1:1:number_of_slices
                        
                        count = count +1;
                        
                        if (E1>1)
                            idx.kspace_encode_step_1(count)=ex.method.PVM_EncSteps1(e1) + round(nY/2);
                        else
                            idx.kspace_encode_step_1(count)=0;  % peu probable
                        end
                        
                        if (E2>1)
                            idx.kspace_encode_step_2(count)=ex.method.PVM_EncSteps2(e2)+ round(nZ/2);
                        else
                            idx.kspace_encode_step_2(count)=0;
                        end                        
                        
                        idx.slice(count)=ns-1;
                        idx.repetition(count)=rep-1;
                        idx.set(count)=nd-1;
                        idx.contrast(count)=ne-1;
                        
                        % idx.phase(count)= ;
                        % idx.set(count)= ;
                        % idx.segment(count)= ;
                        
                        %                 if(
                        %                 idx.parallel_and_calib(count)=
                        %                     if (e1==1)
                        %                         str_msg=sprintf('count %d slice %d e1 %d e2 %d rep %d  ', count, idx.slice(count), idx.kspace_encode_step_1(count), idx.kspace_encode_step_2(count), idx.repetition(count)); disp( str_msg);
                        %                     end
                        
                        % set some flags
                        
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
                        
                        if (E2>1 && e2==1)
                            flag.first_in_encoding_step2(count)=1;
                        else
                            flag.first_in_encoding_step2(count)=0;
                        end
                        
                        if (E2>1 && e2==E2)
                            flag.last_in_encoding_step2(count)=1;
                        else
                            flag.last_in_encoding_step2(count)=0;
                        end
                        
                        if (number_of_repetitions>1 && rep==1)
                            flag.first_in_repetition(count)=1;
                        else
                            flag.first_in_repetition(count)=0;
                        end
                        
                        if (number_of_repetitions>1 && rep==number_of_repetitions)
                            flag.last_in_repetition(count)=1;
                        else
                            flag.last_in_repetition(count)=0;
                        end
                        
                        if (number_of_slices>1 && ns==1)
                            flag.first_in_slice(count)=1;
                        else
                            flag.first_in_slice(count)=0;
                        end
                        
                        if (number_of_slices>1 && ns==number_of_slices)
                            flag.last_in_slice(count)=1;
                        else
                            flag.last_in_slice(count)=0;
                        end
                        
                    end
                end
            end
        end
    end
end

end

