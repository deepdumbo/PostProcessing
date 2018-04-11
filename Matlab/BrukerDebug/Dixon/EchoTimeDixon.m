%% Parameters initialization
    PVM.GradCalCons = 1000.0; % Max Grad Strenght (T/m)

    PVM.Gyro        = 42.577; % Gyromagnetic ratio for water (MHz/T)
    PVM.ChemShift   = 3.57;   % Chemical Shift (ppm)
    PVM.B0          = 400.313032222786; % MHz

    Phase2DGradLim  = 1.0;
    ReadGradLim     = 1.0;
    RewGradLim      = 0.95;

    PVM.EffSWh      = 250000.0; % Bandwidth (Hz)

    PVM.Mat         = [200 200];
    PVM.Fov         = [90.0 90.0];   % Field of view (mm)
    PVM.SpatRes     = PVM.Fov ./ PVM.Mat; % Spatial resolution (mm)

    PVM.MinEchoTime = 2.2; % (ms)
    PVM.EchoTime    = PVM.MinEchoTime;
    PVM.EchoSpacing = 0.0; % Initialize to 0 (ms)

    PVM.AcqTime     = (PVM.Mat(1) / PVM.EffSWh)*1000.0; % Acquisition Time (ms)
    PVM.riseTime    = 0.268; % Gradient rise time (ms)
    PVM.EffPulseDur = 0.9;   % Pulse Duration (ms)
    
    PVM.ReadGradDur = PVM.AcqTime + PVM.riseTime;
    PVM.ReadGrad    = ((1 / (PVM.Fov(1)*0.001)) / ((PVM.Gyro*0.001)*PVM.ReadGradDur)) * 100/PVM.GradCalCons; % ReadGrad in % of max  
    PVM.RewGradDur  = PVM.ReadGradDur;
    PVM.RewGrad     = PVM.ReadGrad;

%% Get First Echo Time
    PVM.StartPt = (1 / (PVM.B0 * PVM.ChemShift)) * 1000.0; % Dixon TE between 2 IP (ms)

    rest    = mod(PVM.EchoTime,PVM.StartPt);
    % Make PVM.EchoTime a multiple of StartPt
    while(rest ~= 0)
        if(mod(PVM.EchoTime - rest, PVM.StartPt) == 0)
            PVM.EchoTime = PVM.EchoTime - rest;
            break;
        elseif(mod(PVM.EchoTime + rest, PVM.StartPt) == 0)
            PVM.EchoTime = PVM.EchoTime + rest;
            break;
        else
            PVM.EchoTime = PVM.EchoTime - rest;
            rest         = mod(PVM.EchoTime,PVM.StartPt);
        end
    end
    
    % Find the shortest PVM.EchoTime matching the minimum requirement
    while(PVM.EchoTime < PVM.MinEchoTime)
        PVM.EchoTime = PVM.EchoTime + PVM.StartPt;
    end
   
%% Update the echo spacing
    % Initialization
    denab           = PVM.riseTime - 0.006;
    minEchoSpacing  = denab + PVM.AcqTime + PVM.riseTime + 0.02 + PVM.RewGradDur;
    
    while(PVM.EchoSpacing < minEchoSpacing)
        PVM.EchoSpacing = PVM.EchoSpacing + (PVM.StartPt / 2.0);
    end
    
    while(mod(PVM.EchoSpacing,PVM.StartPt) == 0)      
        if((PVM.EchoSpacing - (PVM.StartPt / 2.0)) < minEchoSpacing)
            PVM.EchoSpacing = PVM.EchoSpacing + (PVM.StartPt / 2.0);
        else
            PVM.EchoSpacing = PVM.EchoSpacing - (PVM.StartPt / 2.0);
        end
    end
    
    clear RewGradLim ReadGradLim rest denab minEchoSpacing Phase2DGradLim
    