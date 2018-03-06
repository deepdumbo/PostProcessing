% Bloch Equation Simulation
% -----------------------------------------
% 

df = 0;		% Hz off-resonance.
T1 = 650;	% ms.
T2 = 13;	% ms.
TE = 1;		% ms.
dT = 1;
TR = 500;	% ms.
flip = pi/6;	% radians.
Ntr = round(TR/dT);
Nex = 50;	% 20 excitations.

M = [0;0;1];
Rflip = yrot(flip);
[A1,B1] = freeprecess(dT,T1,T2,df);


M(1,Nex*Ntr)=0;	%	Allocate to record all M's.
		% 	Not necessary, but makes program faster.

Mcount=1;
for n=1:Nex
	M(:,Mcount) = Rflip*M(:,Mcount);	

	for k=1:Ntr
		Mcount=Mcount+1;
		M(:,Mcount)=A1*M(:,Mcount-1)+B1;
	end;
end;

time = [0:Mcount-1]*dT;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;