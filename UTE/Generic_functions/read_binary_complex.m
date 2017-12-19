function [ cx_data ] = read_binary_complex( filename, taille )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


fp = fopen(filename);
adata=fread(fp,[taille*2],'float');

for i=1:2:size(adata)
    
    real((i+1)/2,1)=adata(i);
    imag((i+1)/2,1)=adata(i+1);
    
end

for i=1:1:size(real,1)
    
    cx_data(i,1)=complex(real(i,1),imag(i,1));
    
end

% whatever other finalization you need to do
fclose(fp);

end

