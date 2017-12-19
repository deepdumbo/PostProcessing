function [ navigator ] = extract_navigator_resp_and_ecg( data , debut, fin ,repetition , threshold  , Fs_ute )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here





for r=1:1:repetition
    
    %% SELFGATING
    
    % data.navigation.resp=data.filter.resp(debut:fin,:,r);
    % data.navigation.ecg=data.substract.resp(debut:fin,:,r);
    
    debug='N';
    
    offset=0;
    condition=1;
    
    while (condition<=1)
        
        [ matrixC, matrixG, matrixS ] = extract_correlation( data.filter.resp(debut:fin,:,r) , threshold.resp-offset,  debug );
        
        [ matrixU ] = extract_UU( matrixG , debug);
        
        [ navigator.resp.raw(:,r) , coilID] = extract_cluster_clean( matrixU, matrixC, data.filter.resp(:,:,r) , threshold.svd );
        
        condition=sum(coilID);
        
%         str_msg=sprintf('%d %2.3f\n', condition , threshold.resp-offset); disp(str_msg);
        
        offset=offset+0.03;
        
    end
    
    
    offset=0;
    condition=1;
    
    while (condition<=1)
        
        [ matrixC, matrixG, matrixS ] = extract_correlation( data.substract.resp(debut:fin,:,r) , threshold.ecg-offset, debug );
        
        [ matrixU ] = extract_UU( matrixG , debug);
        
        [ navigator.ecg.raw(:,r) , coilID] = extract_cluster_clean( matrixU, matrixC, data.substract.resp(:,:,r) , threshold.svd );
        
        condition=sum(coilID);
        
%         str_msg=sprintf('%d %2.3f\n', condition , threshold.ecg-offset); disp(str_msg);
        
        offset=offset+0.03;
        
    end
    
    %% nous avons extrait les deux navigateurs neanmoins il reste quelques étapes de traitement à effectuer
    
    %% partie respiration
    
    [ navigator.resp.nodrift(:,r) ] = remove_drift_respiration( navigator.resp.raw(:,r), debut, fin );
    
    [ navigator.resp.normalize(:,r) ] = normalize_navigator( navigator.resp.nodrift(:,r), debut, fin );
    
    %% partie ecg
    
    [navigator.ecg.filter(:,r)]=apply_butterworth_filtering(navigator.ecg.raw(:,r),3, Fs_ute);
    
end

end

