function [ input ]=shifting_of_multislice_data_SMS_WS_dim5(input, PE_shift)

% disp('shifting_of_multislice_data_SMS_WS');

center_k_space_sample=round(size(input,2)/2);

number_of_stacks=size(input,4);
slice_acceleration=size(input,5);

if PE_shift~=0
    
for a=1:number_of_stacks 
    for m=1:slice_acceleration
       
     
        
        input(:,:,:,a,m)=(input(:,:,:,a,m)).*...
            repmat(exp(1i*([1:1:size(input,2)]- 1*center_k_space_sample   )*2*pi/PE_shift*(m-1)),[size(input,1) 1  size(input,3) ]  );
%         size(reconSB(:,:,:,nband,:,:,nt))
%        exp( 1i*([1:1:size(reconSB,2)]- 1*center_k_space_sample   )*2*pi/PE_shift*(nband-1));
%         size(exp( 1i*([1:1:size(reconSB,2)]- 1*center_k_space_sample   )*2*pi/PE_shift*(nband-1)));
        
%         size(repmat(exp(1i*([1:1:size(reconSB,2)]- 1*center_k_space_sample   )*2*pi/PE_shift*(nband-1)),[size(reconSB,1) 1 size(reconSB,3) 1 size(reconSB,5) size(reconSB,6)]  ));
    end
end
    
end
          return
