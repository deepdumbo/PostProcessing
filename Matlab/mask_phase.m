m1 = abs(im(:,:,1));

%% 1. Faire un smooth gaussien du 1er echo et le soustraire à l'image original (lisse le bruit)
G = m1 - imgaussfilt(m1);
figure('Name','Echo 1 - smoothed'); imagesc(G);

%% 2. Récupérer std de l'image
st = std(G(:));

% N.B.: l'étape 1 et 2 peuvent être remplacées par la std d'une ROI dans le
% bruit

%% 3. Créer un mask
mask = 0.*m1;
mask(m1 > 3*st) = 1;
figure('Name','Mask initial'); imagesc(mask);

%% 4. Jeu entre erosion / dilatation pour supprimer les pixels parasites et récupérer les pixels manquants
kernel = [0 1 0; 1 1 1; 0 1 0];
% Step 1
er = imerode(mask,kernel);
di = imdilate(er,kernel);

% Step 2
er2 = imerode(di,kernel);

% Step 3
di2 = imdilate(er2,kernel);
er2 = imerode(di2,kernel);


% N.B. : Trouver un moyen de faire ça automatiquement
mask = er2;
figure('Name','Mask ero/dil'); imagesc(mask);


%% 5. Dérouler la phase après avoir appliqué le mask
phi_mask = mask.*(angle(im(:,:,3))-angle(im(:,:,1)));
figure('Name','Phase wrapped mask'); imagesc(phi_mask); colormap gray

unwr = QualityGuidedUnwrap2D_r1(abs(im(:,:,3)).*mask,wrap(phi_mask,pi),0.005);
figure('Name','Phase unwrapped'); imagesc(unwr); colormap gray

line = 140;
figure('Name','Plot comp'); plot(unwr(line,:),'b'); hold on; plot(phi_mask(line,:),'r'); hold on; plot(unwr(line,:)-phi_mask(line,:),'g.-');

legend('Unwrapped phase','Wrapped phase','Unwrapped - Wrapped');

clear line


