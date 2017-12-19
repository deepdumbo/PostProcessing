function [ data_peak_dectection_plus, data_peak_dectection_moins] = find_peaks( input, debut, fin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



data_ecg=input(debut:fin);


%% création de deux listes
data_peak_dectection_plus=0;
data_peak_dectection_moins=0;


for i=2:size(data_ecg)-1;
    
    
    if (data_ecg(i-1)>0 &&  data_ecg(i)<0 )
        
        data_peak_dectection_plus=[data_peak_dectection_plus,i+debut];
        
    end
    
    if (data_ecg(i-1)<0 &&  data_ecg(i)>0 )
        
        data_peak_dectection_moins=[data_peak_dectection_moins,i+debut];
        
    end
end


data_peak_dectection_plus=data_peak_dectection_plus';
data_peak_dectection_plus(:,2)=1;

data_peak_dectection_moins=data_peak_dectection_moins';
data_peak_dectection_moins(:,2)=1;

data_peak_dectection_plus=single(data_peak_dectection_plus);
data_peak_dectection_moins=single(data_peak_dectection_moins);

str_msg=sprintf('Detection de %d pics \n ', size(data_peak_dectection_plus,1)); disp(str_msg);

str_msg=sprintf('Detection de %d pics \n ', size(data_peak_dectection_moins,1));  disp(str_msg);

end

