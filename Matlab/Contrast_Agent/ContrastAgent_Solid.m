
%% Initial parameters
% T1 of the tissu without C.A.
T1tissu     = 1500      / 1000; % s

% T2 of the tissu without C.A.
T2tissu     = 430     / 1000; % s
T2desired   = 50   / 1000; % s

% Contrast Agent
%ca = 'cuso4';
ca = 'mncl2';

if(strcmp(ca,'cuso4'))
    M           = 159.609 * 0.001;  % g/mmol (CuSO4 + 5H2O)
    r1          = 0.67;     % mM/s (CuSO4 + 5H2O)
    r2          = 1.04;     % mM/s (CuSO4 + 5H2O)
end

if(strcmp(ca,'mncl2'))
    M           = 95.211 * 0.001;  % g/mmol (MnCl2 + 5H2O)
    r1          = 6.397;     % mM/s (MnCl2 + 5H2O)
    r2          = 108.266;     % mM/s (MnCl2 + 5H2O)
end

% Volume desired
V           = 350; % mL


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