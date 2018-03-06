
%% Initial parameters
% T1 of the tissu without C.A.
T1tissu = 2.5; % s
% Desired T1
T1desired = 1.05; % s
% Relaxivity rate of the C.A.
r1      = 3.6; % mM/s (Dotarem)
% Initial concentration of the C.A.
Cinit   = 500; % mM -> Dotarem (0.5 mmol/ml)
% Volume of C.A. to be used
V       = 1; % mL



%% Calculation
% Concentration to add to get a specific T1
C = (1/T1desired - 1/T1tissu)/r1; % mM

% Percentage needed
P = (C / Cinit)*100;

% Volume of water needed for the dilution
Vdesired = Cinit*V/C; % mL

fprintf(' Initial T1 : %.2f ms\n Desired T1 : %.2f ms\n Concentration of C.A. needed : %.3f mM (%.3f%% of initial concentration)\n Volume of water for 1x Dilution : %.1f ml (%.1f L)\n',T1tissu*1000, T1desired*1000, C, P, Vdesired, Vdesired / 1000)

if(1)
    T2tissu = 0.230; % s
    r2      = 4.3; % mM/s (Dotarem)
    
    robt = 1/T2tissu + r2*C; % s-1
    T2obt = 1/robt; % s
    
    fprintf('\n It will also affect T2...\n Initial T2 : %.2f ms\n Obtained T2 : %.2f ms\n\n',T2tissu*1000, T2obt*1000)
end