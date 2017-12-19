function [  ] = plot_figure_resp( table, navigator, debut, fin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% 
% close(figure(4))
% figure(4)
% plot(navigator.ecg.testfilter(debut_all:fin_all)); title('ecg navigator'); xlim([1 300000]);
close(figure(6))
figure(6)
subplot(211)


plot(table.respiration.affichage(:,1),'b.');  xlim([debut fin]);   ylim([0 1.5]); 
hold on;
if (size(table.respiration.affichage,2)>1)
plot(table.respiration.affichage(:,2),'c.');
hold on;
end
if (size(table.respiration.affichage,2)>2)
plot(table.respiration.affichage(:,3),'g.');
hold on;
end
if (size(table.respiration.affichage,2)>3)
plot(table.respiration.affichage(:,4),'y.'); 
hold on;
end
if (size(table.respiration.affichage,2)>4)
plot(table.respiration.affichage(:,5),'r.');
hold on;
end
if (size(table.respiration.affichage,2)>5)
plot(table.respiration.affichage(:,6),'r.');
hold on;
end
plot(navigator.allrep.resp.normalize(debut:fin,1)); title('respiratory navigator');
% subplot(212)
% plot(table.respiration_quartile.affichage(:,1),'bo');
% hold on;
% plot(table.respiration_quartile.affichage(:,2),'co');
% hold on;
% plot(table.respiration_quartile.affichage(:,3),'go');
% hold on;
% plot(table.respiration_quartile.affichage(:,4),'yo'); ylim([0 1.5]);
% hold on;
% plot(table.respiration_quartile.affichage(:,5),'ro');
% hold on;
% plot(navigator.allrep.resp.normalize(debut_all:fin_all,1)); title('respiratory navigator'); xlim([debut fin]);


end

