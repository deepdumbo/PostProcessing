%% Works for Bruker rawdata
clear all

%% Initialization
    [ str_user ] = get_PC_name();
    [ str_network_imagerie, str_network_perso ] = get_network_name( str_user );

    % Select the source folder
    patient{1}.folderRaw    = ['/home/',str_user,'/mount/Imagerie/For_Kylian/Dixon/Verification/Ex_Vivo/2D/20180412_Unmasc/' ];
    %patient{2}.folderRaw   = ['/home/',str_user,str_network_imagerie,'DICOM_DATA/2018-03-13_ICM/HIFU_PRECLINIC_18_03_13-11_55_55/' ];
    
    % Select the output folder
    patient{1}.folderOutput = ['/home/',str_user,'/mount/Imagerie/For_Kylian/Dixon/Verification/Ex_Vivo/2D/20180412_Unmasc/'];
    %patient{2}.folderOutput= ['/home/',str_user,str_network_perso,'MatlabUnix/2018-03-13_ARFI/'];
    
    % Select the file name
    patient{1}.filenameOutput = 'README';
    %patient{2}.filenameOutput= 'protocol_arfi_preclinic';

%% Extract info for each patient
    for p=1:size(patient,2)
        
        % Get the current folder
        folderRaw       = patient{p}.folderRaw;
        folderOutput    = patient{p}.folderOutput;
        filenameOutput  = patient{p}.filenameOutput;       

        listing         = dir(folderRaw);
        real_listing    = listing(3:end); % Starts at 3 to avoid '.' and '..' files

        % Display the selected files        
        fprintf(' From : %s \n\t- Size : %d file(s)\n', folderRaw, size(real_listing,1));
        
        % Define the output file
        monFichier      = [folderOutput, filenameOutput, '.txt'];
        monFichierCsv   = [folderOutput, filenameOutput, '.csv'];
        txt_de_sortie   = fopen(monFichier, 'w+');
        csv_de_sortie   = fopen(monFichierCsv, 'w+');
        
        % Fill the output files
        for i=1:size(real_listing,1)
            
            % Pick all the experiment folder paths
            subfolder   = [folderRaw, real_listing(i).name, '/'];          
            
            % Check if the experiment folder is not empty
            if (size(dir(subfolder),1) > 1 && exist([subfolder '/pdata/1/procs'], 'file'))
                
                % Read Bruker files
                info                = read_bru_experiment(subfolder);
                
                % Get the Bruker info
                [ nX, nY, nZ ]              = get_dimensions( info );
                [ nX_acq, nY_acq, nZ_acq ]  = get_dimensions_acq( info );
                [ readout, E1, E2 ]         = get_encoding_size( info , nZ );
                
                % Multi-echo / EchoSpacing managing
                if(info.acqp.NECHOES == 1)
                    EchoSpacing = 'x';
                else
                    EchoSpacing = num2str(info.method.EchoSpacing);
                end
                
                % Headers : Only for the first subfolder
                if (i==1)
                    
                    % Header for the .txt
                    fprintf(txt_de_sortie, '%s \n', info.acqp.ACQ_time{1}(1:10));
                    fprintf(txt_de_sortie, '%s , %s , %5d \n', info.acqp.ORIGIN{1}, info.acqp.ACQ_institution{1} , info.acqp.BF1/42.577  );
                    fprintf(txt_de_sortie, '\n');
                    fprintf(txt_de_sortie, '---------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf(txt_de_sortie, '| %-2s | %-3s | %-12s | %-12s | %-25s | %-3s | %-5s | %-6s | %-2s | %-5s | %-3s | %-10s | %-2s | %-8s |\n' , ' n°','Seq','Start Time','Scan Duration', ' Protocol ','  TR ',' TE  ','   ES  ' ,'ne', ' FA ' , '   Matrix  ','Thickness' , '  Fov  ');
                    fprintf(txt_de_sortie, '---------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                    
                    % Header for the .csv
                    fprintf(csv_de_sortie, ' %-2s , %-3s , %-12s , %-12s , %-25s , %-3s , %-5s , %-6s , %-2s , %-5s , %-3s , %-10s , %-2s , %-8s \n' , ' n°','Seq','Start Time','Scan Duration', ' Protocol ','  TR ',' TE  ','   ES  ' ,'ne', ' FA ' , '   Matrix  ','Thickness' , '  Fov  ' );
                    fprintf(csv_de_sortie, '\n');
                end

                % fprintf(flux_de_sortie, 'n° %d , %42s , %5d , %5d, %5d , %3d \n',i, str_msg, nCoils ,nbSlices  , nReps , nContrasts );
                fprintf(txt_de_sortie, '|  %-2d | %-3d | %-12s | %-13s |  %-24s | %-5.0f | %-5.2f | %-7s | %-2d | %-5d |  %-10s | %-10.2f | %-7s |\n', ...
                    i, str2double(real_listing(i).name), info.acqp.ACQ_time{1}(12:19), info.method.PVM_ScanTimeStr{1}, info.method.Method{1}, info.method.PVM_RepetitionTime, info.method.PVM_EchoTime, EchoSpacing, info.acqp.NECHOES, info.acqp.ACQ_flip_angle, [num2str(nX),'x',num2str(nY),'x',num2str(nZ)], info.method.PVM_SliceThick ,[num2str(info.method.PVM_Fov(1)),'x',num2str(info.method.PVM_Fov(2))]);

                fprintf(csv_de_sortie, '  %-2d , %-3d , %-12s , %-13s ,  %-24s , %-5.0f , %-5.2f , %-7s , %-2d , %-5d ,  %-10s , %-10.2f , %-7s \n', ...
                    i, str2double(real_listing(i).name), info.acqp.ACQ_time{1}(12:19), info.method.PVM_ScanTimeStr{1}, info.method.Method{1}, info.method.PVM_RepetitionTime, info.method.PVM_EchoTime, EchoSpacing, info.acqp.NECHOES, info.acqp.ACQ_flip_angle, [num2str(nX),'x',num2str(nY),'x',num2str(nZ)], info.method.PVM_SliceThick ,[num2str(info.method.PVM_Fov(1)),'x',num2str(info.method.PVM_Fov(2))]);
            
            else
                fprintf('\n Warning : %s is empty\n', subfolder);
            end

        end

        fprintf(txt_de_sortie, '---------------------------------------------------------------------------------------------------------------------------------------------------------\n');

        fclose(txt_de_sortie);
        fclose(csv_de_sortie); 

    end
                 


