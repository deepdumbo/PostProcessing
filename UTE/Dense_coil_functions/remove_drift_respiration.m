function [ output, fitting ] = remove_drift_respiration( navigator_resp, debut, fin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



	 x = linspace(debut, fin, fin-debut+1)';
	  	
	 y=navigator_resp(debut:fin,1);
   		
	 sum_x=0;
	 sum_y=0;
	 sum_xy=0;
	 sum_x2=0;
	
	for  i = 1:1:size(x,1)
	
	 sum_x = sum_x + x(i);
	 sum_y = sum_y+ y(i); 
	 sum_xy = sum_xy+ x(i) * y(i);
	 sum_x2 =sum_x2+  x(i) * x(i); 
    end 

	 mean_x = sum_x / size(x,1);
	 mean_y = sum_y / size(x,1);

	 varx = sum_x2 - sum_x * mean_x;
	 cov = sum_xy - sum_x * mean_y;

	 % check for zero varx

	  param0 = cov / varx;
	 param1 = mean_y - param0 * mean_x;
		
	fitting=zeros(size(navigator_resp));
	
	for  i = 1:1:size(fitting,1)
	
	 fitting(i)=param0*i+param1;
    end
		
	 output=navigator_resp-fitting;
	

end

