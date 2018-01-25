
lala=abs(im(:,:,3));

min_threshold=max(lala(:))*0.05;

fat_threshold=max(lala(:))*0.15;

indice_mask= find( lala> min_threshold );
indice_fat= find( lala> fat_threshold );
% indice_mask_background= find( lala<= min_threshold);

lili=lala;
lulu=lala;
lili(indice_mask(:))=1;
lulu(indice_fat(:))=2;
% lili(indice_mask_background(:))=1;

figure()
subplot(131); imagesc(lala); colormap gray
subplot(132); imagesc(lili); colormap gray
subplot(133); imagesc(lulu); colormap gray
