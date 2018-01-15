S0 = im(:,:,1);
S1 = im(:,:,2);
S2 = im(:,:,3);

close all; 
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
    phi_2 = angle( S2_ );
    %phi_2 = GoldsteinUnwrap2D_r1(abs(S2_),angle(S2_));
   %phi_2 = QualityGuidedUnwrap2D_r1(abs(S2_),angle(S2_));
    phi   = phi_2 / 2;
    
    % Initialize A
    A = sqrt( abs( S2 ) ./ S0_ );

% If we reduce the 3PD to a 2PD we could get W2 = (abs(S0) + abs(S1))/2
% and F2 = (abs(S0) - abs(S1))/2 and abs(S0)+abs(S1)> abs(S0)-abs(S1) then
% W2 will always be brighter than F2 but if the sample has F>W then W2 and
% F2 will be interverted. To avoid this we need to tell when W>F and when
% W<F thanks to the phase as p = +1 for W>F and p = -1 for W<F but for
% better results we will take it as a continuous way.
    pc = cos(angle(S1_.*exp(-1i.*phi)));
    %pc = real(S1_.*exp(-1i.*phi))./abs(S1_);

    % Get W and F images
    W = (S0_ + (pc.*abs(S1))./ A) ./ 2;
    F = (S0_ - (pc.*abs(S1))./ A) ./ 2;
    

% figure, colormap(gray(256))
% imagesc(angle(S2)); 
% xlabel('Pixels'), ylabel('Pixels')
% title('Unwrapped phase image using the 2D-SRNCP algorithm')
% figure 
% surf(double(angle(S2)),'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong')
% view(-30,30), camlight left, axis tight, title('Unwrapped phase image using the 2D-SRNCP displayed as a surface')
% xlabel('Pixels'), ylabel('Pixels'), zlabel('Phase in radians')
     
     % Plot reconstructed images
%      ismrm_imshow(abs(im),[],[1 3],{'echo 1','echo 2','echo 3'}, 'Mag - Wrapped');
%      ismrm_imshow(angle(im),[],[1 3],{'echo 1','echo 2','echo 3'}, 'Phase - Wrapped');
     ismrm_imshow(phi_2,[],[1 1],[], '2phi unwrapped');
     
     % Plot Dixon
     Dix(:,:,1) = W;
     Dix(:,:,2) = F;
     ismrm_imshow(abs(Dix),[],[1 2],{'Water','Fat'}, 'Dixon - Wrapped');