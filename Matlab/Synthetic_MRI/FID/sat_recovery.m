% Bloch Equation Simulation
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
flip = pi/3;	% radians.


M = [0;0;1];
Rflip = yrot(flip);
[Ate,Bte] = freeprecess(TE,T1,T2,df);


M = Rflip*M;	% Magnetization after tip.
M = Ate*M+Bte;	% Magnetization at TE.