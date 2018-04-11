%%% This script is to generate a synthetic FLASH sequence image %%%
% AUTHOR : Kylian HALIOT, PhD - 18/12/2017

disp('--> GENERATE_SYNTHETIC_FLASH_IMAGES');
    disp('    :: Initialization');
    
    close all
    
    addpath('../BrukerDebug/Dixon/');
    
    %% Initialization of the variables
    
    EchoTimeDixon;
    
    % ----------- IRM ------------ %
    Gyro        = 42.577; % Gyromagnetic ratio for water (MHz/T)
    B0          = 400.313032222786; % MHz
    
    % ---------- Matrix ---------- %
    encX   = PVM.Mat(1);
    encY   = PVM.Mat(2);
    encZ   = 1;
    ne     = 5;

    synFLASH = zeros(encX, encY, encZ, ne);

    % --------- Sequence ---------- %
    seqmode = 'FLASH';
    
    TR          = 50;   % ms
    T1_W        = 2500;   % ms
    T1_F        = 500;   % ms
    FA          = (acos(exp(-TR/T1_F))*180)/pi;  % °  
    ChemShift   = 3.57;   % Chemical Shift (ppm)
    dF          = B0*ChemShift;          % Hz
    StartPt     = PVM.StartPt;  % ms

    TE          = PVM.EchoTime;   % ms
    
    if(strcmpi(seqmode,'Dixon'))
        ES          = PVM.EchoSpacing;   % ms
    elseif(strcmpi(seqmode,'FLASH'))
        ES          = StartPt / 2;  % ms
    end
    
    if(ne > 1)
        for i = 2:ne
            TE(i) = TE(i-1) + ES;
        end
    end

    T2e_W  = 130;    % ms
    T2e_F  = 26;    % ms

    S0_W   = 10000;           % scaling factor for water signal
    S0_F   = 6000 ;           % scaling factor for fat signal

%% Initialization of the FLASH signal equations
    disp('    :: Generation of Water and Fat samples');
    for i = 1:ne
        for z = 1:encZ

            % Water fraction
            for y = floor(encY/10):floor(9*encY/10)
                for x = floor(encX/10 - 1):floor(5*encX/10)

                    synFLASH(x,y,z,i) = (S0_W * y * 2) * (((1-exp(-TR / T1_W)) / (1-cos(FA)*exp(-TR/T1_W))) * sin(FA)) * exp(-TE(i)/T2e_W);

                end
            end
    
            % Fat fraction
            for y = floor(encY/10):floor(9*encY/10)
                for x = floor(5*encX/10 + 1):floor(9*encX/10 + 1)

                    synFLASH(x,y,z,i) = (S0_F * y * 2) * (((1-exp(-TR / T1_F)) / (1-cos(FA)*exp(-TR/T1_F))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(-1i*dF*TE(i));

                end
            end
        end
    end

%% Display
TE = TE';
titles = cellstr(strcat('TE = ',num2str(TE(1:size(squeeze(synFLASH),3)),'%-.2f'), ' ms'));
                ismrm_imshow(abs(squeeze(synFLASH)),[],[1 size(squeeze(synFLASH),3)],titles,'Synthetic MRI :: No noise Magnitude');
                ismrm_imshow(angle(squeeze(synFLASH)),[],[1 size(squeeze(synFLASH),3)],titles,'Synthetic MRI :: No noise Phase');  

disp('<-- GENERATE_SYNTHETIC_FLASH_IMAGES');