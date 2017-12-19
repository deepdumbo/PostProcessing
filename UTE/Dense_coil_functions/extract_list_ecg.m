function [ output ] = extract_list_ecg( input, nb_phase_ecg ,repetition_for_reco, nbProjections  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


idxECG=input;

mea=median(diff(idxECG));
stdi=std(diff(idxECG));


%% SORT ECG DATA

 %g.NumberofPhase    %NB of phases
output=zeros(repetition_for_reco* nbProjections, 1);

mea=median(diff(idxECG));
stdi=std(diff(idxECG));
% figure;plot(diff(idxECG))

conteurEx=0;
for i=1:length(output)
    picBefore=idxECG(find(idxECG<i,1,'last'));
    picAfter=idxECG(find(idxECG>=i,1,'first'));
    
    if isempty(picBefore)
        picBefore=idxECG(1)-(idxECG(10)-idxECG(5))/5;
    end
    
    if isempty(picAfter)
        picAfter=idxECG(end)+(idxECG(end)-idxECG(end-3))/3;
    end
    
    idxImg=1+floor((i-picBefore)/(picAfter-picBefore)*nb_phase_ecg);
    
    if idxImg>nb_phase_ecg
        idxImg=nb_phase_ecg;
    elseif idxImg<1
        idxImg=1;
    end
    
    if abs((picAfter-picBefore)-mea)<stdi*3 %& (abs(i-RespBefore)>RespDist) & (abs(i-RespAfter)>RespDist)
        output(i)=idxImg;
    else
        conteurEx=conteurEx+1;
    end
    
end
% exclusions
disp(strcat('ECG data excluded  :', num2str(conteurEx/length(output)*100),'%'));

for i=1:nb_phase_ecg
    lala=find(output==i);
    str_msg=sprintf('%d %d \n',i ,size(lala,1)); disp(str_msg);
end

output=single(output);

end

