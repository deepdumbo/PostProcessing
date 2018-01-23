%%% This script is to generate a synthetic FLASH sequence image %%%
% AUTHOR : Kylian HALIOT, PhD - 18/12/2017

disp('--> GENERATE_SYNTHETIC_FLASH_IMAGES');
    disp('    :: Initialization');
    
    close all
    %% Initialization of the variables
    % ---------- Matrix ---------- %
    encX   = 200;
    encY   = 200;
    encZ   = 1;
    ne     = 30;

    synFLASH = zeros(encX, encY, encZ, ne);

    % --------- Sequence ---------- %
    TR     = 100;   % ms
    T1     = 2500;   % ms
    FA     = (acos(exp(-TR/T1))*180)/pi;  % °
    TE     = 2.8;   % ms
    ES     = 3.15;   % ms

    if(ne > 1)
        for i = 2:ne
            TE(i) = TE(i-1) + ES;
        end
    end

    T2e_W  = 17;    % ms
    T2e_F  = 9;    % ms

    S0_W   = 10000;           % scaling factor for water signal
    S0_F   = 8000 ;           % scaling factor for fat signal

    dF     = 1400;          % Hz

%% Initialization of the FLASH signal equations
    disp('    :: Generation of Water and Fat samples');
    for i = 1:ne
        for z = 1:encZ

            % Water fraction
            for y = floor(encY/4):floor(3*encY/4)
                for x = floor(encX/4):floor(2*encX/4)

                    synFLASH(x,y,z,i) = (S0_W * y * 2) * (((1-exp(-TR / T1)) / (1-cos(FA)*exp(-TR/T1))) * sin(FA)) * exp(-TE(i)/T2e_W);

                end
            end
    
            % Fat fraction
            for y = floor(encY/4):floor(3*encY/4)
                for x = floor(2*encX/4 + 1):floor(3*encX/4)

                    synFLASH(x,y,z,i) = (S0_F * y * 2) * (((1-exp(-TR / T1)) / (1-cos(FA)*exp(-TR/T1))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(-1i*dF*TE(i));

                end
            end
        end
    end

%% Display
addpath('../../ismrm_sunrise_matlab/');

                ismrm_imshow(abs(squeeze(synFLASH)),[],[1 size(squeeze(synFLASH),3)],{});
                ismrm_imshow(angle(squeeze(synFLASH)),[],[1 size(squeeze(synFLASH),3)],{});  

disp('<-- GENERATE_SYNTHETIC_FLASH_IMAGES');