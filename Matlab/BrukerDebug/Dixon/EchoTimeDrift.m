%% Quantifying the impact of little error in the value of Echo Time for Dixon
%% In phase and Out of phase state

    B0          = 400.313032222786; % MHz
    ChemShift   = 3.57; % ppm
    dF          = B0 * ChemShift; % Hz
    StartPt     = (1 / dF) * 1000.0; % Dixon TE between 2 IP (ms)
    EchoTime    = 0; % (ms)

    nTour       = 1; % Number of 2pi cycles (odd value expected)
    nEchoes     = 100; % Number of desired echoes

    
%% Echo Time Drift
    EchoSpacing = (StartPt / 2) * nTour; % (ms)
    rEchoSpacing= ceil(EchoSpacing*100)/100; % EchoSpacing used in Bruker
    
    AllEchoTimes(1)  = 0;
    rAllEchoTimes(1) = 0;
    for i=2:nEchoes
        AllEchoTimes(i) = AllEchoTimes(i-1) + EchoSpacing;
        rAllEchoTimes(i) = rAllEchoTimes(i-1) + rEchoSpacing;
    end
    
    diff = abs(rAllEchoTimes - AllEchoTimes);
    for i=1:nEchoes   
        Integrated(i) = sum(diff(1:i));
    end
    
%% Phase Drift
    AllPhase  = 2*pi*dF*(AllEchoTimes/1000);
    rAllPhase = 2*pi*dF*(rAllEchoTimes/1000);
    
    pdiff = abs(rAllPhase - AllPhase);
    for i=1:nEchoes   
        Dephasing(i) = sum(pdiff(1:i));
    end
    
    Drift = [AllEchoTimes; rAllEchoTimes; diff; Integrated; AllPhase; rAllPhase]';
    
 %% Plot
 figure('Name','EchoTimeDrift.m :: Error Plot');
    subplot(121);
        plot(1:nEchoes, Integrated, 'b*'); 
        title('Cumulative Error at each Echo Time');
        xlabel('Echo number');
        ylabel('Cumulative Error (in ms)');
    subplot(122);
        plot(rAllEchoTimes, Dephasing/pi,'b*');
        title('Dephasing Error at each Echo Time (used in Bruker)');
        xlabel('Echo Time (in ms)');
        ylabel('Cumulative Dephasing (in number of \pi)');