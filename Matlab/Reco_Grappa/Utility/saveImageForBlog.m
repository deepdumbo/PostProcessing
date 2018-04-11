function [ output ] = saveImageForBlog( image_complexe ,folder, name, ext)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nCoils=size(image_complexe,3);



for c = 1:nCoils
    magnitude(:,:,c)=abs(image_complexe(:,:,c));
    phase(:,:,c)=angle(image_complexe(:,:,c));   
    
end

rotation=0;

if (strcmp(ext,'.dat'))
    
    for c=1:1:nCoils
        str_c=sprintf('%d',c);
        
        filename=[folder,'magnitude_', name, str_c,ext]
        matrix=imrotate(magnitude(:,:,c),rotation);
        save(filename,'matrix','-ascii');
        
        filename=[folder,'phase_', name, str_c,ext ]
        matrix=imrotate(phase(:,:,c),rotation);
        save(filename,'matrix','-ascii');
        
    end;
end;


% if (ext=='.png')
%     
%     for s=1:1:dimz
%         str_s=sprintf('%d',s);
%         filename=[folder, name, str_s,ext ]
%         matrix=imrotate(input(:,:,s,1),rotation);
%         imwrite(matrix,filename);
%         
%     end;
%     
% end;




output=1;
end

