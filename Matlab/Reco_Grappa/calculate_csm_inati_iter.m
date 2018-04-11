function [ coil_map_out, coil_combined_out ] = calculate_csm_inati_iter( img_tempo, smoothing_value, niter, thresh )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%    def calculate_csm_inati_iter(im, smoothing=5, niter=5, thresh=1e-3,
%                              verbose=False):
%     """ Fast, iterative coil map estimation for 2D or 3D acquisitions.
%     Parameters
%     ----------
%     im : ndarray
%         Input images, [coil, y, x] or [coil, z, y, x].
%     smoothing : int or ndarray-like
%         Smoothing block size(s) for the spatial axes.
%     niter : int
%         Maximal number of iterations to run.
%     thresh : float
%         Threshold on the relative coil map change required for early
%         termination of iterations.  If ``thresh=0``, the threshold check
%         will be skipped and all ``niter`` iterations will be performed.
%     verbose : bool
%         If true, progress information will be printed out at each iteration.
%     Returns
%     -------
%     coil_map : ndarray
%         Relative coil sensitivity maps, [coil, y, x] or [coil, z, y, x].
%     coil_combined : ndarray
%         The coil combined image volume, [y, x] or [z, y, x].
%     Notes
%     -----
%     The implementation corresponds to the algorithm described in [1]_ and is a
%     port of Gadgetron's ``coil_map_3d_Inati_Iter`` routine.
%     For non-isotropic voxels it may be desirable to use non-uniform smoothing
%     kernel sizes, so a length 3 array of smoothings is also supported.
%     References
%     ----------
%     .. [1] S Inati, MS Hansen, P Kellman.  A Fast Optimal Method for Coil
%         Sensitivity Estimation and Adaptive Coil Combination for Complex
%         Images.  In: ISMRM proceedings; Milan, Italy; 2014; p. 4407.
%

verbose=1

img_tempo2=permute(img_tempo , [3,2,1]);

%     im = np.asarray(im)

if (ndims(img_tempo2) < 3 || ndims(img_tempo2) > 4 )
    disp('error');
end

if (ndims(img_tempo2) == 3)
    % pad to size 1 on z for 2D + coils case
    images_are_2D = 1;
    im(:, 1, :, :) = img_tempo2;
else
    images_are_2D = 0;
end

size(im)

%     convert smoothing kernel to array
%     if isinstance(smoothing, int):

smoothing = zeros( 3,1)
smoothing(:) = smoothing_value

ndims(smoothing)

if (ndims(smoothing) > 1 || size(smoothing,1) ~= 3)
    disp('smoothing should be an int or a 3-element 1D array');
end

%     if (images_are_2D==1)
%         smoothing[2] = 1  # no smoothing along z in 2D case
%
%     % smoothing kernel is size 1 on the coil axis
%     smoothing = np.concatenate(([1, ], smoothing), axis=0)

ncha = size(im,1);

%     try:
%         % numpy >= 1.7 required for this notation
%         D_sum = im.sum(axis=(1, 2, 3))
D_sum=sum(sum(sum(im,4),3),2);

%     except:
%         D_sum = im.reshape(ncha, -1).sum(axis=1)

%     v = 1/np.linalg.norm(D_sum)
v=1/ norm(D_sum);
%     D_sum *= v
D_sum=D_sum*v;
R = 0;

for cha =1:1:ncha
    %         R += np.conj(D_sum[cha]) * im[cha, ...]
    R= R + conj(D_sum(cha) * im(cha, :,:,:));
end



%     eps = np.finfo(im.real.dtype).eps * np.abs(im).mean()
for it =1:1:niter
    if verbose
        disp(sprintf('Coil map estimation: iteration %d of %d', it, niter));
    end
    
    if (thresh > 0)
        %             prevR = R.copy()
        prevR = R;
    end
    
    %         R = np.conj(R)
    R = conj(R);
    %         coil_map = im * R[np.newaxis, ...]
    for cha =1:1:ncha
        RR(cha,1,:,:)=R(1,1,:,:);
    end
    coil_map = im .* RR;
    %         coil_map_conv = smooth(coil_map, box=smoothing)
    h_smooth=ones(smoothing_value)/(smoothing_value^2);
    
    % Smooth the covariance
    for cha =1:1:ncha
        coil_map_conv(cha,1,:,:) = smooth_function(squeeze(coil_map(cha,1,:,:)), h_smooth ,smoothing);
    end
    %         D = coil_map_conv * np.conj(coil_map_conv)
    D = coil_map_conv .* conj(coil_map_conv);
    %         R = D.sum(axis=0)
    R = sum(D,1);
    %         R = np.sqrt(R) + eps
    R = sqrt(R) + eps;
    %         R = 1/R
    R = 1/R;
    %         coil_map = coil_map_conv * R[np.newaxis, ...]
    for cha =1:1:ncha
        RR(cha,1,:,:)=R(1,1,:,:);
    end
    coil_map=coil_map_conv.*RR;
    %         D = im * np.conj(coil_map)
    D = im .* conj(coil_map);
    %         R = D.sum(axis=0)
    R = sum(D,1);
    %         D = coil_map * R[np.newaxis, ...]
    for cha =1:1:ncha
        RR(cha,1,:,:)=R(1,1,:,:);
    end
    D = coil_map .* RR;
    
    %         try:
    %             # numpy >= 1.7 required for this notation
    %             D_sum = D.sum(axis=(1, 2, 3))
    D_sum=sum(sum(sum(D,4),3),2);
    %         except:
    %             D_sum = im.reshape(ncha, -1).sum(axis=1)
    
    %         v = 1/np.linalg.norm(D_sum)
    v = 1/norm(D_sum);
    
    D_sum = D_sum*v;
    %
    imT = 0;
    for cha=1:1:ncha
        %             imT += np.conj(D_sum[cha]) * coil_map[cha, ...]
        imT = imT+ conj(D_sum(cha) * coil_map(cha,:,:,:));
    end
    %             magT = np.abs(imT) + eps
    magT = abs(imT) + eps;
    %         imT /= magT
    imT = imT./ magT;
    %         R = R * imT
    R = R .* imT;
    %         imT = np.conj(imT)
    imT = conj(imT);
    %         coil_map = coil_map * imT[np.newaxis, ...]
    for cha =1:1:ncha
        imTT(cha,1,:,:)=imT(1,1,:,:);
    end
    
    coil_map = coil_map .* imTT;
    
    if (thresh > 0)
       
        diffR = squeeze(R) - squeeze(prevR);    
        vRatio = norm(diffR,2) / norm(squeeze(R),2);
        if (verbose==1)
            disp(vRatio)
        end
        if (vRatio < thresh)
            break
        end
        
    end
    
        coil_combined = squeeze(sum((im .* conj(coil_map)),1));
    
    
end

if (images_are_2D==1)
    % remove singleton z dimension that was added for the 2D case
    %         coil_combined = coil_combined[0, :, :]
    %         coil_map = coil_map[:, 0, :, :]
end

size(coil_combined)
size(coil_map)

coil_combined_out=permute(coil_combined,[2 1]);
coil_map_out=permute(squeeze(coil_map), [3 2 1]);

size(coil_combined_out)
size(coil_map_out)
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