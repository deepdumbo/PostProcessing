function [ W, F, IP, OP ] = Dixon_3P( S0, S1, S2, e )
% This function aims to give a water/fat-only images with 3 points (echos)
% with a phase encoding of the chemical shift as following (0, pi, 2pi).
% Parameters :
%   S0  : first echo at theta = 0
%   S1  : second echo at theta = pi
%   S2  : third echo at theta = 2pi
%   W   : Water-only image
%   F   : Fat-only image
%   IP  : In-Phase image
%   OP  : Out-of-Phase image
%   e   : Amplitude loss error tolerated

% BASED ON : Multipoint Dixon Technique for Water and Fat Proton and 
% Susceptibility Imaging, Glover G. - J Magn Reson Imaging 1991;1:521?530.

% AUTHOR : Kylian HALIOT, PhD - 16/10/2017

%% Check arguments
    narginchk( 3, nargin('Dixon_3P') );
    if nargin == 3
        e = 0;
    end
    if e > 1
        e = 1;
    end

%% In-Phase and Out-of-Phase images + Amplitude loss factor
    IP = zeros(size(S0));
    OP = zeros(size(S0));
    theta0 = angle( S0 );

    S0_ = abs( S0 );
    S1_ = S1 .* exp( -1i.*theta0 );
    S2_ = S2 .* exp( -1i.*theta0 );

    Two_theta = angle( S2_ ) - angle( S0_ ); %afficher
    theta = Two_theta ./ 2;
    A_square = abs( S2 ) ./ S0_ ;
    A = sqrt( A_square ); % Amplitude loss factor

    % Calculate IP and OP
    IP( A >= 1 - e ) = ( S0_(A >= 1-e) + abs( S2(A >= 1-e) ) ) ./ 2;
    IP( A < 1 - e ) = S0_(A < 1-e) ;

    OP( A >= 1 - e ) = S1_(A >= 1-e) .* exp(-1i.*theta(A >= 1-e));
    OP( A < 1 - e ) = ( S1_(A < 1-e) .* exp(-1i.*theta(A < 1-e)) ) ./ A(A < 1-e);

%% Water/Fat-only images reconstruction
    W = (IP + OP) ./ 2;
    F = (IP - OP) ./ 2;


end

