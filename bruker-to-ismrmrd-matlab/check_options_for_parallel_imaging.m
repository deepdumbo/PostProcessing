function [ Ysamp_regular , Ysamp_ACS , acceleration_factor_y ] = check_options_for_parallel_imaging( ex )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


%% get some info with readable name

[ nX, nY, nZ ] = get_dimensions( ex );

[ nX_acq, nY_acq, nZ_acq ] = get_dimensions_acq( ex );


if (nY~=nY_acq  && ex.method.PVM_EncPft(2)==1 )
    
    acceleration_factor_y=2;
    ACShw = 12 ;
    nombre_de_ligne_de_acs_x=ex.method.PVM_EncPpiRefLines(1);
    nombre_de_ligne_de_acs_y=ex.method.PVM_EncPpiRefLines(2);
    fprintf('nombre_de_ligne_de_acs %d \n', nombre_de_ligne_de_acs_y);  
    fprintf('acceleration suivant y  %d facteur \n', acceleration_factor_y); 
    
else
    
   acceleration_factor_y=1;
   ACShw = 0 ;
   fprintf('pas acceleration suivant y  facteur %d \n', acceleration_factor_y); 
    
end


lala= ex.method.PVM_EncSteps1+ round(nY/2);
lili=lala(2:end)-lala(1:end-1);

u=find(lili==1);
v=find(lili==2);

% size(lala)

% lala(u+1)
% lala(v)

% set ACS lines for GRAPPA simulation (fully sampled central k-space
% region)
% GRAPPA ACS half width i.e. here 28 lines are ACS


%% TODO a generaliser  ( le facteur 2 )

if (acceleration_factor_y>1)
    
    Ysamp_u = [1:2:nY] ; % undersampling by every alternate line
    Ysamp_regular=Ysamp_u-1;
    Ysamp_ACS = [nY/2-ACShw+1 : nY/2+ACShw] ; % GRAPPA autocalibration lines
    Ysamp = union(Ysamp_u, Ysamp_ACS)-1 ; % actually sampled lines
    nYsamp = length(Ysamp) ; % number of actually sampled
    
    for  i = 1:1:length(Ysamp)
        
%         str_smg=sprintf('%d %d ', Ysamp(i),  ex.method.PVM_EncSteps1(i) + round(nY/2)); disp(str_smg);
        
        if (Ysamp(i)==ex.method.PVM_EncSteps1(i) + round(nY/2))
            
        else
            fprintf('%d %d \n', Ysamp(i),  ex.method.PVM_EncSteps1(i) + round(nY/2)); %disp(str_smg);
            disp('pb encodage y a a cause de l imagerie parallele');
        end
        
    end
    
    disp('fin gestion parallel imaging');
    
    size(u);
    size(v);
    size(u)+size(v);
    lala(u);
    Ysamp_ACS;
    
else
    
    Ysamp_regular=[];
    Ysamp_ACS=[];
    
end


end

