function [ dest ] = ButterworthLowpassFilterSixthOrder( src , a, b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NZEROS = 6;
NPOLES = 6;

GAIN =  a(1)/b(1);

xv=zeros(NZEROS+1,1) ;
yv=zeros(NPOLES+1,1);

for i = 1:1:size(src,1)
    
    %xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; xv[4] = xv[5]; xv[5] = xv[6];
    xv(1) = xv(2); xv(2) = xv(3); xv(3) = xv(4); xv(4) = xv(5); xv(5) = xv(6); xv(6) = xv(7);
    %xv[6] = src[i] / GAIN;
    xv(7) = src(i,1) / GAIN;
    %yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; yv[4] = yv[5]; yv[5] = yv[6];
    yv(1) = yv(2); yv(2) = yv(3); yv(3) = yv(4); yv(4) = yv(5); yv(5) = yv(6); yv(6) = yv(7);
    
    %         yv[6] =   (xv[0] + xv[6]) + 6.0 * (xv[1] + xv[5]) + 15.0 * (xv[2] + xv[4])
    %                      + 20.0 * xv[3]
    %                      + ( -0.1396600417 * yv[0]) + (  1.1086708553 * yv[1])
    %                      + ( -3.7230194289 * yv[2]) + (  6.7850160254 * yv[3])
    %                      + ( -7.0995038188 * yv[4]) + (  4.0616439992 * yv[5]);
    %         dest[i] = yv[6];
    
    %         yv(7) =   (xv(1) + xv(7)) + 6.0 * (xv(2) + xv(6)) + 15.0 * (xv(3) + xv(5))  + 20.0 * xv(4)  ...
    %                      + ( -0.1396600417 * yv(1)) + (  1.1086708553 * yv(2)) ...
    %                      + ( -3.7230194289 * yv(3)) + (  6.7850160254 * yv(4)) ...
    %                      + ( -7.0995038188 * yv(5)) + (  4.0616439992 * yv(6));
    
    
    yv(7) =   (xv(1) + xv(7)) + 6.0 * (xv(2) + xv(6)) + 15.0 * (xv(3) + xv(5))  + 20.0 * xv(4)  ...
        + ( -a(7) * yv(1)) + ( -a(6) * yv(2))  ...
        + ( -a(5) * yv(3)) + ( -a(4) * yv(4)) ...
        + (-a(3) * yv(5)) + ( -a(2) * yv(6));
    
    
    dest(i,1) = yv(7);
end



end


