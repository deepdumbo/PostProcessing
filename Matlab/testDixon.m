S0 = im(:,:,1);
S1 = im(:,:,2);
S2 = im(:,:,3);

close all;

%% T2* mapping (requires at least 10-15 echoes)
disp('--> T2* Mapping');
    % Get the T2* map of the dataset
    ld = 0;
   
    if(ld == 0)
        disp('    :: Fat T2* mapping');
        Fout = T2star_mapping_bruker(im,TE, 'Fat area');
        disp('    :: Water T2* mapping');
        Wout = T2star_mapping_bruker(im,TE, 'Water area');
        %save('T2stars.mat', 'Fout', 'Wout');
    else
        tmp = load('T2star.mat');
        Fout = tmp.Fout;
        Wout = tmp.Wout;
    end
    
    % As the mapping is done depending on a user-defined ROI, take only the
    % mean of the T2* map
%     Fout.fitmap(Fout.fitmap == 0) = NaN;
%     fT2star                       = mean(Fout.fitmap(:), 'omitnan');
%     Wout.fitmap(Wout.fitmap == 0) = NaN;
%     wT2star                       = mean(Wout.fitmap(:), 'omitnan');
      fT2star                       = Fout.fitmean;
      wT2star                       = Wout.fitmean;
    
    % Plot of the fit for water and fat
    pixel = 1;
    TEc = linspace(0,TE(end),500); % Continuous TE for better plotting the fit
    Ffat = Fout.res(pixel,1) + Fout.res(pixel,2).*exp(-TEc./Fout.res(pixel,3));
    Fwat = Wout.res(pixel,1) + Wout.res(pixel,2).*exp(-TEc./Wout.res(pixel,3));

disp('<-- T2* Mapping');  

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
disp('--> Dixon_3P with T2*');

%% 0. Keep only the first 3 TEs
    disp('    :: Parameters initialization');
    TE_  = TE(1:3);
    dTE = TE_(2) - TE_(1);

%% Other way to find W and F by saying R = -1/T2*
%% 1. Remove phi0 from the equations
    % Get Rw and Rf
    Rw = -1 / wT2star;
    Rf = -1 / fT2star;
    
    % Get phi0
    phi0 = angle( S0 );

    % Remove phi0 from the equations and Rw*t0 which is constant in all
    % equations
    S0_ = abs( S0 ) .* exp(-Rw*TE_(1));
    S1_ = S1 .* exp( -1i.*phi0 ).* exp(-Rw*TE_(1)).*exp(-Rw*dTE);
    S2_ = S2 .* exp( -1i.*phi0 ).* exp(-Rw*TE_(1)).*exp(-2*Rw*dTE);
    
%% 2. Get phi and calculate W and F images
% We're using 2phi from S2 instead of phi from S1 because S1 is in an
% opposed phase state for W and F and the SNR is lower than S2 which is
% in an in phase state for W and F. 
    
    % Get phi
    disp('    :: Phase masking & unwrapping');
%      tresh = 0.005; % Treshold for the noise
%      phi_2 = QualityGuidedUnwrap2D_r1(abs(S2_),wrap(angle(S2_),pi), tresh);
     phi_2 = Masked_Unwrap(S0, S2_, 1);
     phi   = phi_2 / 2;
    
    % Calculate A (which is a ratio between the 2 T2* factors)
    A2 = exp(2*(Rf-Rw)*dTE);
    A  = sqrt(A2);

    % If we reduce the 3PD to a 2PD we could get W2 = (abs(S0) + abs(S1))/2
    % and F2 = (abs(S0) - abs(S1))/2 and abs(S0)+abs(S1)> abs(S0)-abs(S1) then
    % W2 will always be brighter than F2 but if the sample has F>W then W2 and
    % F2 will be interverted. To avoid this we need to tell when W>F and when
    % W<F thanks to the phase as p = +1 for W>F and p = -1 for W<F but for
    % better results we will take it as a continuous way.
    disp('    :: Water & Fat attribution map');
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

%% 3. Get rid of phi and Rw*dTE
    S2__ = abs(S2_);
    S1__ = abs(S1_);

%% 4. Get W and F images
    disp('    :: Water & Fat separation');
    W = (S2__ + (A - A2).*S0_ + pc.*S1__) ./ (2 + A - A2);
    F = (S2__ - pc.*S1__) ./ ((A2+A).*exp((Rf-Rw).*TE_(1)));

disp('<-- Dixon_3P with T2*');
%% SNR calculation
     clear Dix shw
     Dix(:,:,1) = W;
     Dix(:,:,2) = F;
    
%     SNR.im1     = SNR_ROI('Noise Area', im(:,:,1), Wout.fitmask);
%     SNR.W       = SNR_ROI('',Dixon(:,:,1), Wout.fitmask, SNR.im1.noise.mask);
%     SNR.F       = SNR_ROI('',Dixon(:,:,2), Fout.fitmask, SNR.im1.noise.mask);
%     SNR.Wstar   = SNR_ROI('',Dix(:,:,1), Wout.fitmask, SNR.im1.noise.mask);
%     SNR.Fstar   = SNR_ROI('',Dix(:,:,2), Fout.fitmask, SNR.im1.noise.mask);
    

    
%% Plots
disp('--> Figure plots');
     % Plot T2* fitting
     disp('    :: T2* fitted curves');
        figure('Name','Water and Fat : T_{2}^{*} fitting','NumberTitle','off');
        subplot(211); 
            plot(TE,Fout.SFull(:,pixel),'o'); 
            hold on; 
            plot(TEc,Ffat); 
            title(['Fat (T_{2}^{*} = ',num2str(Fout.fitmean),' ms)']);   
            legend('Experimental points','F = x_0 + S\ite^{-t / T_{2}^{*}}');
            xlabel('t (ms)');
            ylabel('Signal magnitude');     
        subplot(212); 
            plot(TE,Wout.SFull(:,pixel),'o'); 
            hold on; 
            plot(TEc,Fwat); 
            title(['Water (T_{2}^{*} = ',num2str(Wout.fitmean),' ms)']); 
            legend('Experimental points','F = x_0 + S\ite^{-t / T_{2}^{*}}');
            xlabel('t (ms)');
            ylabel('Signal magnitude');
            
    % Plot unwrapped phase
    disp('    :: Unwrapped phase');
    clear phases
    phases(:,:,1) = angle(S2_);
    phases(:,:,2) = phi_2;
    ismrm_imshow(phases,[],[1 2],{'2\phi wrapped' '2\phi unwrapped'},'Phase unwrapping');
    
    % Plot hc and pc 
    disp('    :: Water & Fat attribution maps');
    tt(:,:,1) = hc;
    tt(:,:,2) = pc;
    fig_unwr  = ismrm_imshow(tt,[],[1 2],{'hc' 'pc'},'Phase unwrapping');
                     
    % Plot all methods (Dixon without and with T2* + diff
    disp('    :: Water & Fat images');
     shw(:,:,1) = Dixon(:,:,1);
     shw(:,:,2) = Dix(:,:,1);
     shw(:,:,3) = Dixon(:,:,2);
     shw(:,:,4) = Dix(:,:,2);
    
%     titles = {['Water no T_{2}^{*} (SNR = ',num2str(SNR.W.val),')']  ...
%               ['Water T_{2}^{*} (SNR = ',num2str(SNR.Wstar.val),')'] ...
%               ['Fat no T_{2}^{*} (SNR = ',num2str(SNR.F.val),')']    ...
%               ['Fat T_{2}^{*} (SNR = ',num2str(SNR.Fstar.val),')']};  
    titles = {'Water no T_{2}^{*}'  ...
              'Water T_{2}^{*}'     ...
              'Fat no T_{2}^{*}'    ...
              'Fat T_{2}^{*}'}; 
    fig_cmp = ismrm_imshow(abs(shw),[],[2 2],titles,'Comparison');
 
disp('<-- Figure plots');