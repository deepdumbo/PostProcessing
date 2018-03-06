
%% Initial parameters
% T1 of the tissu without C.A.
T1tissu     = 1.43; % s
T1desired   = 0.5; % s

% T2 of the tissu without C.A.
T2tissu     = 0.230; % s
T2desired   = 0.1; % s

% Contrast Agent
M           = 159.609 * 0.001;  % g/mmol (CuSO4 + 5H2O)
r1          = 0.67;     % mM/s (CuSO4 + 5H2O)
r2          = 1.04;     % mM/s (CuSO4 + 5H2O)

% Volume desired
V           = 50; % mL


%% Calculation based on T2
% Concentration to add to get a specific T2
C = (1/T2desired - 1/T2tissu)/r2; % mM

% Corresponding mass
m = M * C * V;

fprintf(' Initial T2 : %.2f ms\n Desired T2 : %.2f ms\n Mass of C.A. needed : %.3f mg (%.3f mM)\n',T2tissu*1000, T2desired*1000, m, C)

if(1)
  
    robt = 1/T1tissu + r1*C; % s-1
    T1obt = 1/robt; % s
    
    fprintf('\n It will also affect T1...\n Initial T1 : %.2f ms\n Obtained T1 : %.2f ms\n\n',T1tissu*1000, T1obt*1000)
end