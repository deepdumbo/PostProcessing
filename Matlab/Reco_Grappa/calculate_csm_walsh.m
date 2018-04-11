function [ csm_out, rho ] = calculate_csm_walsh( img_tempo, smoothing, niter )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     def calculate_csm_walsh(img, smoothing=5, niter=3):
%     '''Calculates the coil sensitivities for 2D data using an iterative version of the Walsh method
%     :param img: Input images, ``[coil, y, x]``
%     :param smoothing: Smoothing block size (default ``5``)
%     :parma niter: Number of iterations for the eigenvector power method (default ``3``)
%     :returns csm: Relative coil sensitivity maps, ``[coil, y, x]``
%     :returns rho: Total power in the estimated coils maps, ``[y, x]``
%     '''

img=permute(img_tempo , [3,2,1]);

size(img)

%     assert img.ndim == 3, %"Coil sensitivity map must have exactly 3 dimensions"

ncoils = size(img,1);
ny = size(img,2);
nx = size(img,3);

% Compute the sample covariance pointwise
Rs = zeros(ncoils,ncoils,ny,nx);


for p=1:1:ncoils
    for q=1:1:ncoils
        Rs(p,q,:,:) = img(p,:,:) .* conj(img(q,:,:));
    end
end

h_smooth=ones(smoothing)/(smoothing^2);

% Smooth the covariance
for p=1:1:ncoils
    for q=1:1:ncoils
        Rs(p,q,:,:) = smooth_function(squeeze(Rs(p,q,:,:)), h_smooth ,smoothing);
    end
end

% At each point in the image, find the dominant eigenvector
% and corresponding eigenvalue of the signal covariance
% matrix using the power method
rho = zeros(ny, nx);
csm = zeros(ncoils, ny, nx);
for y=1:1:ny
    for x=1:1:nx
        %             R = Rs[:,:,y,x]
        R = Rs(:,:,y,x);
        %             v = np.sum(R,axis=0)
        v = sum(R,1);
        %             lam = np.linalg.norm(v)
        lam=norm(v);  %lam=squeeze(sqrt(sum(real(v).^2 + imag(v).^2,2)));
        v = v'/lam;
        
        %
        for iter =1:1:niter            
            v = R*v;
            lam=norm(v);   % norm(v)
            v = v/lam;
        end
        rho(y,x) = lam;
        csm(:,y,x) = v;
    end
end

csm_out=permute(csm , [3,2,1]);


return


function [ smooth_img ] = smooth_function(img, h_smooth, box)
%     '''Smooths coil images
%     :param img: Input complex images, ``[y, x] or [z, y, x]``
%     :param box: Smoothing block size (default ``5``)
%     :returns simg: Smoothed complex image ``[y,x] or [z,y,x]``
%     '''

smooth_img=conv2(img,h_smooth,'same');

%     t_real = zeros(size(img))
%     t_imag = zeros(size(img))
%
%     ndimage.filters.uniform_filter(img.real,size=box,output=t_real)
%     ndimage.filters.uniform_filter(img.imag,size=box,output=t_imag)
%
%     simg = t_real + 1j*t_imag

return
