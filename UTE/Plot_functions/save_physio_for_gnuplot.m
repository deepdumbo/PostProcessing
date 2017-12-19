function [  ] = save_physio_for_gnuplot( input_1, input_2 , folder, filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

datatosave(:,1)=input_1;
datatosave(:,2)=input_2;

st_msg=sprintf('min %f max %f ', min(datatosave(:,2)) , max(datatosave(:,2))); disp(st_msg);
disp([folder,filename ]);
save([folder,filename ],'datatosave','-ascii');

end

