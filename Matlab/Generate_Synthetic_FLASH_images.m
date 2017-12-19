%%% This script is to generate a synthetic FLASH sequence image %%%
% AUTHOR : Kylian HALIOT, PhD - 18/12/2017

%% Initialization of the variables
% ---------- Matrix ---------- %
encX   = 200;
encY   = 200;
encZ   = 1;
ne     = 3;

synFLASH = zeros(encX, encY, encZ, ne);

% --------- Sequence ---------- %
TR     = 100  / 1000;   % ms
T1     = 2500 / 1000;   % ms
FA     = acos(-TR/T1);  % °
TE     = 2.8  / 1000;   % ms
ES     = 3.15  / 1000;   % ms

if(ne > 1)
    for i = 2:ne
        TE(i) = TE(i-1) + ES;
    end
end

T2e_W  = 50  / 1000;    % ms
T2e_F  = 110 / 1000;    % ms

S0_W   = 100;           % scaling factor for water signal
S0_F   = 50 ;           % scaling factor for fat signal

dF     = 1400;          % Hz

%% Initialization of the FLASH signal equations

for i = 1:ne
    for z = 1:encZ
        
        % Water fraction
        for y = floor(encY/4):floor(3*encY/4)
            for x = floor(encX/5):floor(2*encX/5)
                
                synFLASH(x,y,z,i) = (S0_W * y * 2) * (((1-exp(-TR / T1)) / (1-cos(FA)*exp(-TR/T1))) * sin(FA)) * exp(-TE(i)/T2e_W);
                
            end
        end
        
        % Fat fraction
        for y = floor(encY/4):floor(3*encY/4)
            for x = floor(3*encX/5):floor(4*encX/5)
                
                synFLASH(x,y,z,i) = (S0_F * y * 2) * (((1-exp(-TR / T1)) / (1-cos(FA)*exp(-TR/T1))) * sin(FA)) * exp(-TE(i)/T2e_F) * exp(1i*2*pi*dF*TE(i));
                
            end
        end
    end
end

%% Display
for s=1:ne
    subplot(2,ne,s);        imagesc(abs(squeeze(synFLASH(:,:,s))));     axis square;    title(['Echo', num2str(s)]);  
    subplot(2,ne,s+ne);     imagesc(angle(squeeze(synFLASH(:,:,s))));   axis square;
    colormap(gray);
end   
