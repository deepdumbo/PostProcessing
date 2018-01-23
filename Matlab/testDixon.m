S0 = im(:,:,1);
S1 = im(:,:,2);
S2 = im(:,:,3);

close all;

%% T2* mapping (requires at least 10-15 echoes)
    % Get the T2* map of the dataset
    
    ld = 0;
    
    if(ld == 0)
        Fout = T2star_mapping_bruker(im,TE, 'Fat area');
        Wout = T2star_mapping_bruker(im,TE, 'Water area');
        %save('T2stars.mat', 'Fout', 'Wout');
    else
        tmp = load('T2star.mat');
        Fout = tmp.Fout;
        Wout = tmp.Wout;
    end
    
    % As the mapping is done depending on a user-defined ROI, take only the
    % mean of the T2* map
    Fout.fitmap(Fout.fitmap == 0) = NaN;
    fT2star                       = mean(Fout.fitmap(:), 'omitnan');
    Wout.fitmap(Wout.fitmap == 0) = NaN;
    wT2star                       = mean(Wout.fitmap(:), 'omitnan');
    
    % Plot of the fit for water and fat
    pixel = 1;
    Ffat = Fout.res(pixel,1) + Fout.res(pixel,2).*exp(-Fout.TE./Fout.res(pixel,3));
    Fwat = Wout.res(pixel,1) + Wout.res(pixel,2).*exp(-Wout.TE./Wout.res(pixel,3));
    figure('Name','Water and Fat : T_{2}^{*} fitting','NumberTitle','off');
        subplot(211); 
            plot(Fout.TE,Fout.S(:,pixel),'o'); 
            hold on; 
            plot(Fout.TE,Ffat); 
            title(['Fat (T_{2}^{*} = ',num2str(Fout.res(pixel,3)),' ms)']);   
            legend('Experimental points','F = x_0 + S\ite^{-TE / T_{2}^{*}}');
            xlabel('TE (ms)');
            ylabel('Signal magnitude');     
        subplot(212); 
            plot(Wout.TE,Wout.S(:,pixel),'o'); 
            hold on; plot(Wout.TE,Fwat); 
            title(['Water (T_{2}^{*} = ',num2str(Wout.res(pixel,3)),' ms)']); 
            legend('Experimental points','F = x_0 + S\ite^{-TE / T_{2}^{*}}');
            xlabel('TE (ms)');
            ylabel('Signal magnitude');
    
%% We have basically 6 equations (3 magnitude and 3 phase)
%  - phi = phase accumulated because of the B0 inhomogeneity during the
%          echo shift delta.
%  - phi0 = all other phase errors (gradient, receive chain group delays, 
%           eddy currents, concomitant fields, receive field B1 phase and
%           for GRE also the phase accumulated from B0 inhomogeneity
%           during the TE of the in-phase scan.
%  - Aw Af= amplitude loss for 3 sources (diffusion effect, T2* effect and 
%           spectral broadening) but Glover and al. 1991 considers
%           diffusion effect as negligble
%   Equations to resolve are then :
%   - S0 = (W.*Aw0 + F.*Af0).*exp(i*phi0)
%   - S1 = (W.*Aw1 - F.*Af1).*exp(i(phi0 + phi))
%   - S2 = (W.*Aw2 + F.*Af2).*exp(i(phi0 + 2phi))

%% 0. Keep only the first 3 TEs
    TE_  = TE(1:3);
    dTE = TE_(2) - TE_(1);
    
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
     %phi_2 = angle( S2_ );
    phi_2 = QualityGuidedUnwrap2D_r1(abs(S2_),angle(S2_),0.01);
    phi   = phi_2 / 2;
    
    clear phases
    phases(:,:,1) = angle(S2_);
    phases(:,:,2) = phi_2;
    ismrm_imshow(phases,[],[1 2],{'2\phi wrapped' '2\phi unwrapped'},'Phase unwrapping');
    
    % Initialize Tw and Tf
    Tw0   = exp(((wT2star - fT2star) ./ (wT2star*fT2star)) .* TE_(1));
    dTw   = exp(((wT2star - fT2star) ./ (wT2star*fT2star)) .* dTE  );
    
    Tf0   = exp(((fT2star - wT2star) ./ (wT2star*fT2star)) * TE_(1));
    dTf   = exp(((fT2star - wT2star) ./ (wT2star*fT2star)) * dTE  );   

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

clear x y

%     % Get W and F images according to a fatT2star factorisation
    W = (((pc.*abs(S1)-abs(S2)).*exp(dTE/fT2star) + S0_.*(1+ exp(-dTE/fT2star))) ./  (Tf0.*(1 + dTf + (1 - dTf.^2).*exp(-dTE/fT2star)))) .* exp(TE_(1)./fT2star);
            
    F = (((pc.*abs(S1)-abs(S2)).*exp(dTE/fT2star) - S0_.*dTf.*(1 - dTf).*exp(-dTE/fT2star)) ./ (1 + dTf + (1 - dTf.^2).* exp(-dTE/fT2star)) ) .* exp(TE_(1)./fT2star);

      % Get W and F images according to a fatT2star factorisation
%       W = (((pc.*abs(S1)-abs(S2)).*exp(dTE/wT2star) + S0_.*dTw.*(1 + dTw.*exp(-dTE/wT2star))) ./ (1 + dTw + (dTw.^2 - 1).* exp(-dTE/wT2star)) ) .* exp(TE_(1)./wT2star);
%       
%       F = ((S0_.*(1 - exp(-dTE/wT2star)) - (pc.*abs(S1)-abs(S2)).*exp(dTE/wT2star)) ./  (Tw0.*(2 + (dTw - 1).*exp(-dTE/wT2star)))) .* exp(TE_(1)./wT2star);

%% SNR calculation
     clear Dix diff shw
     Dix(:,:,1) = W;
     Dix(:,:,2) = F;
     diff = abs(Dixon) - abs(Dix);
    
    SNR.im1     = SNR_ROI('Noise Area', im(:,:,1), Wout.fitmask);
    SNR.W       = SNR_ROI('',Dixon(:,:,1), Wout.fitmask, SNR.im1.noise.mask);
    SNR.F       = SNR_ROI('',Dixon(:,:,2), Fout.fitmask, SNR.im1.noise.mask);
    SNR.Wstar   = SNR_ROI('',Dix(:,:,1), Wout.fitmask, SNR.im1.noise.mask);
    SNR.Fstar   = SNR_ROI('',Dix(:,:,2), Fout.fitmask, SNR.im1.noise.mask);
    
     
     % Plot all methods (Dixon without and with T2* + diff)
     
     shw(:,:,1) = Dixon(:,:,1);
     shw(:,:,2) = Dix(:,:,1);
     shw(:,:,3) = diff(:,:,1);
     shw(:,:,4) = Dixon(:,:,2);
     shw(:,:,5) = Dix(:,:,2);
     shw(:,:,6) = diff(:,:,2);
    
    titles = {['Water no T_{2}^{*} (SNR = ',num2str(SNR.W.val),')']  ...
              ['Water T_{2}^{*} (SNR = ',num2str(SNR.Wstar.val),')'] ...
               'Diff Water'                                          ...
              ['Fat no T_{2}^{*} (SNR = ',num2str(SNR.F.val),')']    ...
              ['Fat T_{2}^{*} (SNR = ',num2str(SNR.Fstar.val),')']   ...
               'Diff Fat' };                                    
    ismrm_imshow(abs(shw),[],[2 3],titles,'Comparison')