function [ W, F ] = Dixon_3P( S0, S1, S2, TE)
% This function aims to give a water/fat-only images with 3 points (echos)
% with a phase encoding of the chemical shift as following (0, pi, 2pi).
% Parameters :
%   S0  : first echo at theta = 0
%   S1  : second echo at theta = pi
%   S2  : third echo at theta = 2pi
%   W   : Water-only image
%   F   : Fat-only image
%   TE  : Echo times

% BASED ON : Multipoint Dixon Technique for Water and Fat Proton and 
% Susceptibility Imaging, Glover G. - J Magn Reson Imaging 1991;1:521?530.

% AUTHOR : Kylian HALIOT, PhD - 16/10/2017

%% Check arguments
    narginchk( 4, nargin('Dixon_3P') ); 
    
%% We have basically 6 equations (3 magnitude and 3 phase) with a 7th that 
%  constraints the others. And 5 unknowns (W, F, phi0, phi and A):
%  - phi = phase accumulated because of the B0 inhomogeneity during the
%          echo shift delta.
%  - phi0 = all other phase errors (gradient, receive chain group delays, 
%           eddy currents, concomitant fields, receive field B1 phase and
%           for GRE also the phase accumulated from B0 inhomogeneity
%           during the TE of the in-phase scan.
%  - A    = amplitude loss for 3 sources (diffusion effect, T2* effect and 
%           spectral broadening) but Glover and al. 1991 considers
%           diffusion effect as negligible and because fat isn't composed
%           by only 1 chemical specie its T2 is a mean of all the species
%           and shorter than water's then it has been decided to focus on
%           the T2' component of the T2*.
%   Equations to resolve are then :
%   - S0 = (W+F).*exp(i*phi0)
%   - S1 = (W-F).*A.*exp(i(phi0 + phi))
%   - S2 = (W+F).*A².*exp(i(phi0 + 2phi))
%   - with angle(S0) + angle(S2) = 2*angle(S1)


%% 1. Remove phi0 from the equations
    % Get phi0
    phi0 = angle( S0 );

    % Remove phi0 from the equations
    S0_ = abs( S0 );
    S1_ = S1 .* exp( -1i.*phi0 );
    S2_ = S2 .* exp( -1i.*phi0 );
    
%% 2. Get phi and calculate W and F images
% We're using 2phi from S2 instead of phi from S1 because S1 is in an
% opposed phase state for W and F and the SNR is lower than S2 which is
% in an in phase state for W and F. 
    
    % Get phi
    % phi_2 = angle( S2_ );
     phi_2 = QualityGuidedUnwrap2D_r1(abs(S2_),angle(S2_), 0.05);
     phi   = phi_2 / 2;
    
    ismrm_imshow(phi_2,[],[1 1],[], '2phi unwrapped');
    
    % Initialize A
    A = sqrt( abs( S2 ) ./ S0_ );

% If we reduce the 3PD to a 2PD we could get W2 = (abs(S0) + abs(S1))/2
% and F2 = (abs(S0) - abs(S1))/2 and abs(S0)+abs(S1)> abs(S0)-abs(S1) then
% W2 will always be brighter than F2 but if the sample has F>W then W2 and
% F2 will be interverted. To avoid this we need to tell when W>F and when
% W<F thanks to the phase as p = +1 for W>F and p = -1 for W<F but for
% better results we will take it as a continuous way.
hc = cos(angle(S1_.*exp(-1i.*phi)));
pc = NaN(size(hc));

for x = 1:size(hc, 2)
    for y = 1:size(hc,1)
        
        if(0.5 <= hc(x,y) && hc(x,y) <= 1)
            pc(x,y) = 1;
        elseif ( -0.5 < hc(x,y) && hc(x,y) < 0.5)
            pc(x,y) = 0;
        elseif (-1 <= hc(x,y) && hc(x,y) <= -0.5)
            pc(x,y) = -1;
        end
        
    end
end
   
    % Get W and F images
    W = (S0_ + (pc.*abs(S1))./ A) ./ 2;
    F = (S0_ - (pc.*abs(S1))./ A) ./ 2;


end

