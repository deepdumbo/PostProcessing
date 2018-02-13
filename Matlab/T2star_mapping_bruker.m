function out = T2star_mapping_bruker(im, TE, name)

    %Prog fitting T2star mappping from GRE Bruker
    %T2star_mapping_bruker.m

    %% Initialization
        % Don't pick the first echo which is not in steady state yet
        %TE      = TE(2:end)';
        TE      = TE(1:2:end)';
        I       = abs(im(:,:,1:2:end)); % Use for the fit on In-phase state
        Ifull   = abs(im); % Use for recovery all the echoes and for plotting against the fit curve
        image   = I(:,:,1);

    %% Création du masque - Mise en forme des data
        % Show the image of interest
        RoiFig = figure('Name','Please draw a ROI...' , 'NumberTitle','off'); 
        imagesc(abs(im(:,:,1))); title(upper(name)); colormap(gray); colorbar;   zoom off;

        % Draw the ROI 
        fitmask_signal = double(roipoly);

        % Convert the data
        numb_pixels_analyzed_signal=sum(sum(fitmask_signal));

        % Close the fig
        close(RoiFig);

        %S: vector of all pixels in the ROI during 1 TE on 2 (Fit only the
        %In-phase echoes)
        tmp = 0 ;
        S   = [];
        for k = 1:1:size(I,3)
            for i = 1:1:size(I,1)
                for j = 1:1:size(I,2)
                        if (fitmask_signal(i,j) == 1) 
                            tmp     = tmp + 1 ;
                            S(k,tmp)= I(i,j,k); 
                        end
                end
            end
            tmp = 0;
        end
        
        %Full: vector of all pixels in the ROI during TE
        tmp        = 0 ;
        FullData   = [];
        for k = 1:1:size(Ifull,3)
            for i = 1:1:size(Ifull,1)
                for j = 1:1:size(Ifull,2)
                        if (fitmask_signal(i,j) == 1) 
                            tmp     = tmp + 1 ;
                            FullData(k,tmp)= Ifull(i,j,k); 
                        end
                end
            end
            tmp = 0;
        end

    %% Fit
        % Waitbar initialisation
        h = waitbar(0, 'Work in progress...', 'Name','T2* Mapping');

        % Initialisation
        twoD_output_matrix = zeros(numb_pixels_analyzed_signal, 3) ;  
        image_a         = 0*image;
        image_T2star    = 0*image;
        image_b         = 0*image;
        options=optimset('Display','off', 'Algorithm','levenberg-marquardt');

        % Fit
        for k=1:1:numb_pixels_analyzed_signal

            % Update waitbar
            step = k / numb_pixels_analyzed_signal;
            waitbar(step, h , ['Work in progress... (',num2str(k),'/',num2str(numb_pixels_analyzed_signal),')']);

            % Pixel per pixel fitting
            x0(1) = 0;
            x0(2) = S(1,k);
            x0(3) = 0;
            res   = lsqcurvefit(@fitT2, x0, TE(:), S(:,k), [], [], options);

            % Calcul du T2*
            T2star = -1./res(3);
            twoD_output_matrix(k, 1:end) = [res(1), res(2), T2star];

            tp=0;
            for i=1:1:size(image,1)
                for j=1:1:size(image,2)
                    if (fitmask_signal(i,j)==1),
                        tp=tp+1;
                        image_b(i,j)         = twoD_output_matrix( tp,   1  ) ;
                        image_a(i,j)         = twoD_output_matrix( tp,   2  ) ;
                        image_T2star(i,j)    = twoD_output_matrix( tp,   3  ) ;
                    end
                end
            end

               image_T2star(image_T2star<0) = 0;
        end
        close(h);
        
        out.fitmap  = image_T2star;
        out.S       = S;
        out.SFull   = FullData;
        out.res     = twoD_output_matrix;
        out.TE      = TE;
        out.fitmask = fitmask_signal;
end

function F = fitT2(x,grille_echo)
    F = x(1)+x(2)*exp(grille_echo*x(3));
end