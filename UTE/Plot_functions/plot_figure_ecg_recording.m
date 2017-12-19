function [  ] = plot_figure_ecg_recording( t_ute, ecg_resample, locs , T_ute, debut, fin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


close(figure(6))
figure(6)
subplot(1,1,1);
plot(t_ute,-ecg_resample(:,2),'r'); xlim([debut fin]);
hold on;
plot(locs*T_ute, 200, '*');
end

