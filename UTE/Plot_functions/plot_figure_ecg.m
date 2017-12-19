function [  ] = plot_figure_ecg( t_ecg , ecg_truth, t_ute , navigator, data_peak_dectection,T_ute, debut, fin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


close(figure(7))
figure(7)
subplot(1,1,1);
plot(t_ecg,-ecg_truth(:,2)/2-400,'r'); xlim([debut fin]);
hold on;
plot(t_ecg,ecg_truth(:,3)/2+200,'r');  
hold on;
plot(t_ecg,-ecg_truth(:,4)/10,'r'); 
hold on;
plot(t_ute,navigator.ecg.testfilter*2e5-600,'b');
hold on;
plot(data_peak_dectection.plus_new(:,1)*T_ute,data_peak_dectection.plus_new(:,2), '*b');
hold on;
plot(data_peak_dectection.moins_new(:,1)*T_ute,data_peak_dectection.moins_new(:,2), '*g');
hold on;


end

