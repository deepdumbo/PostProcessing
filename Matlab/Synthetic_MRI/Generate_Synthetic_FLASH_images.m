%%% This script is to generate a synthetic FLASH sequence image %%%
% AUTHOR : Kylian HALIOT, PhD - 18/12/2017

disp('--> GENERATE_SYNTHETIC_FLASH_IMAGES');
    disp('    :: Initialization');
    
    close all
    
    addpath('../BrukerDebug/Dixon/');
    
    %% Initialization of the variables
    DEBUG = 1; % 0 : none
               % 1 : Dixon
               % 2 : Dixon + mapping + SNR
    
    EchoTimeDixon;
         
    % ---------- Matrix ---------- %
    encX   = PVM.Mat(1);
    encY   = PVM.Mat(2);
    encZ   = 1;
    ne     = 10;
    synFLASH = zeros(encX, encY, encZ, ne);
    
    c = [encX/2 encY/2];
    r = 30;
    
    % ----------- IRM ------------ %
    Gyro        = PVM.Gyro;                     % (MHz/T)
    B0          = PVM.B0;                       % MHz
    Precession  = 1;                            %  1 : Clockwise
                                                % -1 : Counter-Clockwise

    % --------- Sequence ---------- %
    seqmode = 'flash';
    
    TR          = 50;                            % ms
    T1_W        = 2500;                          % ms
    T1_F        = 500;                           % ms
    T2e_W       = 25;                            % ms
    T2e_F       = 130;                           % ms
    S0_W        = 500 ;                          % scaling factor for water signal
    S0_F        = 250 ;                          % scaling factor for fat signal
    
    FA          = (acos(exp(-TR/T1_F)));         % Ernst angle for fat (rad)  
    ChemShift   = PVM.ChemShift;                 % Chemical Shift (ppm)
    dF          = B0*ChemShift;                  % Hz
    StartPt     = PVM.StartPt;                   % ms
%     TE          = PVM.EchoTime;                % ms
    TE          = 0;
    
    if(strcmpi(seqmode,'Dixon'))
        ES          = PVM.EchoSpacing;           % ms
    elseif(strcmpi(seqmode,'FLASH'))
        pi_mod      = 1;
        ES          = pi_mod * (StartPt / 2);    % ms
    end
    
    if(ne > 1)
        for i = 2:ne
            TE(i) = TE(i-1) + ES;
        end
    end

%% Initialization of the FLASH signal equations
    disp('    :: Generation of Water and Fat samples');
    for i = 1:ne
        for z = 1:encZ

            % Square shaped phantom
            for y = floor(encY*0.1):floor(encY*0.9)
                
                % Water fraction
                for x = floor(encX*0.1 - 1):floor(encX*0.5)
                    synFLASH(x,y,z,i) = (S0_W * y) * (((1-exp(-TR / T1_W)) / (1-cos(FA)*exp(-TR/T1_W))) * sin(FA)) * exp(-TE(i)/T2e_W) * exp(Precession*1i*2*pi*B0*TE(i));
                end 
                
                % Fat fraction
                for x = floor(encX*0.5 + 1):floor(encX*0.9 + 1)
                    synFLASH(x,y,z,i) = (S0_F * y) * (((1-exp(-TR / T1_F)) / (1-cos(FA)*exp(-TR/T1_F))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(Precession*1i*(2*pi*B0 + dF)*TE(i));
                end
            end
            
            % Add local fat areas in water fraction
            for x = 70:80
                for y = 65:75
                    synFLASH(x,y,z,i) = (S0_F * y) * (((1-exp(-TR / T1_F)) / (1-cos(FA)*exp(-TR/T1_F))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(Precession*1i*(2*pi*B0 + dF)*TE(i));
                end
                for y = 125:135
                    synFLASH(x,y,z,i) = (S0_F * y) * (((1-exp(-TR / T1_F)) / (1-cos(FA)*exp(-TR/T1_F))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(Precession*1i*(2*pi*B0 + dF)*TE(i));
                end
            end
            
            for x = 50:60
                for y = 95:105
                    synFLASH(x,y,z,i) = (S0_F * y) * (((1-exp(-TR / T1_F)) / (1-cos(FA)*exp(-TR/T1_F))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(Precession*1i*(2*pi*B0 + dF)*TE(i));
                end
            end
            
            for x = 90:91
                for y = 50:150
                    synFLASH(x,y,z,i) = (S0_F * y) * (((1-exp(-TR / T1_F)) / (1-cos(FA)*exp(-TR/T1_F))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(Precession*1i*(2*pi*B0 + dF)*TE(i));
                end
            end
            
            % Make the phantom round
            for y = 1:encY 
                for x = 1:encX                   
                    if(sqrt((x-c(1)).^2 + (y-c(2)).^2) >= (encX/2 - r))
                        synFLASH(x,y,z,i) = 0;
                    end                   
                end
            end
            
            % Add an extra local fat area not connected to the phantom
            for x = 15:35
                for y = 15:35
                    synFLASH(x,y,z,i) = (S0_F * y) * (((1-exp(-TR / T1_F)) / (1-cos(FA)*exp(-TR/T1_F))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(Precession*1i*(2*pi*B0 + dF)*TE(i));
                end
            end
            
        end
    end

    
%% Debug and mappings
    % Check IP OOP IP
    if(DEBUG == 1)
        for i=1:size(synFLASH,4)
            fprintf(' Echo %d\n   W : %f rad.\n   F : %f rad.\n   Diff : %f rad. (%f x pi)\n',i,angle(synFLASH(100,100,1,i)),angle(synFLASH(101,100,1,i)),(angle(synFLASH(100,100,1,i))-angle(synFLASH(101,100,1,i))),(angle(synFLASH(100,100,1,i))-angle(synFLASH(101,100,1,i)))/pi);
        end
    end
    
    % T2* mapping
    if(DEBUG == 2)
        disp('    :: T2* mapping');
        disp('        :: Fat T2* mapping');
        Fout = T2star_mapping_bruker(im,TE, 'Fat area');
        disp('        :: Water T2* mapping');
        Wout = T2star_mapping_bruker(im,TE, 'Water area');
        
        % Plot of the fit for water and fat
        pixel = 1;
        TEc = linspace(0,TE(end),500); % Continuous TE for better plotting the fit
        Ffat = Fout.res(pixel,1) + Fout.res(pixel,2).*exp(-TEc./Fout.res(pixel,3));
        Fwat = Wout.res(pixel,1) + Wout.res(pixel,2).*exp(-TEc./Wout.res(pixel,3));
        
        disp('        :: T2* fitted curves');
        figure('Name','Water and Fat : T_{2}^{*} fitting','NumberTitle','off');
        subplot(211); 
            plot(TE,Fout.SFull(:,pixel),'o'); 
            hold on; 
            plot(TEc,Ffat); 
            title(['Fat (T_{2}^{*} = ',num2str(Fout.fitmean),' ms)']);   
            legend('Experimental points','F = x_0 + S\ite^{-t / T_{2}^{*}}');
            xlabel('t (ms)');
            ylabel('Signal magnitude');     
        subplot(212); 
            plot(TE,Wout.SFull(:,pixel),'o'); 
            hold on; 
            plot(TEc,Fwat); 
            title(['Water (T_{2}^{*} = ',num2str(Wout.fitmean),' ms)']); 
            legend('Experimental points','F = x_0 + S\ite^{-t / T_{2}^{*}}');
            xlabel('t (ms)');
            ylabel('Signal magnitude');
    end
        
        
%% Display
    TE = TE';
    titles = cellstr(strcat('TE = ',num2str(TE(1:size(squeeze(synFLASH),3)),'%-.2f'), ' ms'));
                    ismrm_imshow(abs(squeeze(synFLASH)),[],[1 size(squeeze(synFLASH),3)],titles,'Synthetic MRI :: No noise Magnitude');
                    ismrm_imshow(angle(squeeze(synFLASH)),[],[1 size(squeeze(synFLASH),3)],titles,'Synthetic MRI :: No noise Phase');  

disp('<-- GENERATE_SYNTHETIC_FLASH_IMAGES');