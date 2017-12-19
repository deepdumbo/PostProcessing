function [  ] = plot_figure_ecg_sort_phase( t_ecg, ecg_truth, t_ute , navigator , data_peak_dectection, T_ute, table_ecg_recording, table_ecg_self_gating_in_systole ,  debut, fin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




close(figure(7))
figure(7)
subplot(1,1,1);
plot(t_ecg,-ecg_truth(:,2)/2-400,'r');
hold on;
plot(t_ute,navigator.ecg.testfilter*2e5-600,'b');  xlim([debut fin])
hold on;
plot(data_peak_dectection.moins_new(:,1)*T_ute,data_peak_dectection.moins_new(:,2)-300, '*g');
hold on;
plot(t_ute, table_ecg_recording*5-300,'r');
hold on;
plot(t_ute, table_ecg_self_gating_in_systole*5-300,'b');

end

