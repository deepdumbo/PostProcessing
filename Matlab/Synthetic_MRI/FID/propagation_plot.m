% Bloch Equation Simulation, Excercise B-1e
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
TR = 500;	% ms.
flip = pi/3;	% radians.


M = [0;0;1];
Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR,T1,T2,df);
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atetr,Btetr] = freeprecess(TR-TE,T1,T2,df);

%	Calculation using B-1d first:

Mss = inv(eye(3)-Atr*Rflip)*Btr;
Mte1 = Ate*Rflip*Mss+Bte


% 	Direct calculation at TE

% 	Starting at TE, M=M1
%	At TR, M=M2, and M2=Atetr*M1+Btetr.
%	At TE, M=M3, and M3=Ate*Rflip*M2+Bte.
%			M3=Ate*Rflip*(Atetr*M1+Btetr)+Bte.
%
%	But M3=M1=Mte2 in steady state:

Mte2 = inv(eye(3)-Ate*Rflip*Atetr)* (Ate*Rflip*Btetr+Bte)