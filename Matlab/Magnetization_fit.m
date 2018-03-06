%% This function aims to calculte the original magnetization of the concerned tissue
function [out] = Magnetization_fit(im,TR,T1w,T1f,FA,TE)
% TR    : Repetition Time in ms
% T1w   : Relaxation time of tissue in ms
% T1f   : Relaxation time of fat in ms
% FA    : Flip Angle in ° 
% TE    : Echo Time in ms
    
    addpath('Synthetic_MRI/');
    
    %% Initialization
    FA = FA * pi/180; % degrees to rad
    I = abs(im(:,:,1)); 
    B0 = 400e6; % Hz
    
    Fout = T2star_mapping_bruker(im,TE,'Fat Area');
    Wout = T2star_mapping_bruker(im,TE,'Water Area');
    out.SNR = SNR_ROI('Noise Area',I);
    
    % Plot of the fit for water and fat
    pixel = 1;
    TEc = linspace(0,TE(end),500); % Continuous TE for better plotting the fit
    Ffat = Fout.res(pixel,1) + Fout.res(pixel,2).*exp(-TEc./Fout.res(pixel,3));
    Fwat = Wout.res(pixel,1) + Wout.res(pixel,2).*exp(-TEc./Wout.res(pixel,3));
    
    % Magnetisation calculation
    i = 0;
    ampli = 1000;
    
    for k = TEc
        Fmag(i+1)  = steadystate_signal(FA,T1f,Fout.fitmean,k,TR,B0,'yes') * ampli;
        Wmag(i+1)  = steadystate_signal(FA,T1w,Wout.fitmean,k,TR,B0,'yes') * ampli;
        i = i+1;
    end
    
    %% Plot
    figure('Name','Water and Fat : T_{2}^{*} fitting','NumberTitle','off');
       subplot(211);
            % Fat
            plot(TE,Fout.SFull(:,pixel),'bo'); 
            hold on; 
            plot(TEc,Ffat,'b'); 
            
            % Water
            hold on; 
            plot(TE,Wout.SFull(:,pixel),'ro'); 
            hold on; 
            plot(TEc,Fwat,'r');
            hold off;
            
            title(['Water (T_{2}^{*} = ',num2str(Wout.fitmean),' ms) & Fat (T_{2}^{*} = ',num2str(Fout.fitmean),' ms)']);  
            legend('Fat area Exp. points','Fat Fitting Curve', 'Water area Exp. points', 'Water Fitting Curve');
            xlabel('t (ms)');
            ylabel('Signal magnitude');
            
        subplot(212);
            % Fat
            plot(TEc,abs(Fmag),'b');
            hold on;

            % Water
            plot(TEc,abs(Wmag),'r');
            hold off;
            
            title(['Magnetization']);  
            legend('Fat area Mag.','Water area Mag.');
            xlabel('t (ms)');
            ylabel('Signal magnitude');

end