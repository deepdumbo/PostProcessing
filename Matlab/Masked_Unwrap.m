function [unwr, mask] = Masked_Unwrap(S0, S2_, Debug)

    %% 1. Faire un smooth gaussien du 1er echo et le soustraire à l'image original (lisse le bruit)
    m1 = abs(S0);
     
%     %% 2. Récupérer std de l'image
%     [noiseROI,~] = ROImask(m1,'Noise Area');
%     st = std(noiseROI(:),'omitnan');
%     G = m1 - st;
%     
%     %% 3. Créer un mask
%     mask = 0.*m1;
%     mask(m1 > 3*st) = 1;
% 
%     %% 4. Jeu entre erosion / dilatation pour supprimer les pixels parasites et récupérer les pixels manquants
%     kernel = [0 1 0; 1 1 1; 0 1 0];
%     % Step 1
%     er = imerode(mask,kernel);
%     di = imdilate(er,kernel);
% 
%     % Step 2
%     er2 = imerode(di,kernel);
% 
%     % Step 3
%     di2 = imdilate(er2,kernel);
%     er2 = imerode(di2,kernel);
%     
%     % N.B. : Trouver un moyen de faire ça automatiquement
%     mask = er2;
      tresh = 0.1;
      idx = m1 > tresh*max(m1(:));
      mask = 0.*m1;
      mask(idx) = 1;

    %% 5. Dérouler la phase après avoir appliqué le mask
    % tresh = 0.005; %default tresh ! Used before
    angS2_ = angle(S2_);
    unwr = QualityGuidedUnwrap2D_r1(abs(S2_).*mask,angS2_.*mask,tresh);
    %unwr = phaseUnwrap2D(angS2_.*mask);
    
    if(Debug)
        figure('Name','Echo 1 - smoothed'); imagesc(G);
        figure('Name','Mask initial'); imagesc(mask);
        figure('Name','Mask ero/dil'); imagesc(mask);
        figure('Name','Phase unwrapped'); imagesc(unwr); colormap gray
        line = 70;
        figure('Name','Plot comp'); plot(unwr(line,:),'b'); hold on; plot(angS2_(line,:),'r'); hold on; plot(wrap(unwr(line,:)-angS2_(line,:),0),'g.-');
        legend('Unwrapped phase','Wrapped phase','Unwrapped - Wrapped');
        clear line
    end
end