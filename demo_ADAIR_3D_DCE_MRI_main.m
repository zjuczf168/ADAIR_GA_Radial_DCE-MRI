%--------------------------------------------------------------------------
%% load data, Cartesian has been tested, non-Cart should be normalized, R2017
%--------------------------------------------------------------------------
clear all;close all; %#ok<CLALL>

set(0,'DefaultFigureWindowStyle','docked')

addpath('.\grasp_v2\');
addpath(genpath('.\Matlab_Toolbox_czf\BM3D_MRI_toolbox\'))

addpath('.\grasp_v2\nufft_toolbox\');
% define number of spokes to be used per frame (Fibonacci number)
nspokes=21;%%spoke can be any Fibonacci No.
% load DCE-MRI data
load breast_data.mat
b1=b1/max(abs(b1(:)));
% data dimensions
[nx,ntviews,nc]=size(kdata);
for ch=1:nc,kdata(:,:,ch)=kdata(:,:,ch).*sqrt(w);end
% number of frames
nt=floor(ntviews/nspokes);
% 
kdata=kdata(:,1:nt*nspokes,:);
k=k(:,1:nt*nspokes); 
w=w(:,1:nt*nspokes);
% sort the data into a time-series
for ii=1:nt
    kdatau(:,:,:,ii)=kdata(:,(ii-1)*nspokes+1:ii*nspokes,:);
    ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
    wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
end
% multicoil NUFFT operator, can also be used for Cartesian
param.E=MCNUFFT(ku,wu,b1);
% undersampled data
param.y=kdatau;
clear kdata kdatau k ku wu w
% nufft recon
recon_nufft=param.E'*param.y; 
% parameters for reconstruction             
param.W = TV_Temp();
lambda = 0.02*max(abs(recon_nufft(:)));%%lambda

img_slc = param.E'*param.y;%ifft2call(kspace_slc);

% img_r1 = coil_combine(img_slc, sens_slc, 3);

[N(1), N(2), num_chan] = size(img_slc);


%--------------------------------------------------------------------------
%% SENSE
%--------------------------------------------------------------------------

% apply subsampling mask
kspace_raw = param.y;%kspace_slc .* mask_coils;

img_zf = flipdim(recon_nufft,1);%coil_combine( ifft2call(kspace_coils), sens_slc ); 

num_iter = 20;  


img_pocs = img_zf;

for t = 1:num_iter

    tmp_k = param.E*img_pocs;
    img_pocs = img_pocs + param.E'*(kspace_raw-tmp_k);

end

% mosaic(img_pocs, 1, 1, fig_num, ['pocs sense: ', num2str(rmse(img_pocs, img_r1))], [0,.4])%, setGcf(.5)

  
%--------------------------------------------------------------------------
%% denoising-MRI
%--------------------------------------------------------------------------

num = 30;                   % outer iters

inner_iter_num1 = 1;       
inner_iter_num2 = 10;

inner_iter_num = logspace(log10(inner_iter_num1), log10(inner_iter_num2), num); 
 
sigma_bm3ds = 36;      % initial 
final_noise = 1;        % final 

sigma_bm3d = logspace(log10(sigma_bm3ds),log10(final_noise),num);


I2 = double(kspace_raw);

% vBM3D-MRI iterations
tic
I11 = img_zf;

for kp = 1:num  
    disp(num2str(sigma_bm3d(kp)))        
    
        for inner = 1:inner_iter_num(kp) 
%             
            min_r = 1*min(real(I11(:)));
            I11_r = real(I11) - min_r;
            scl_r = max(I11_r(:));
            I11_r = I11_r /scl_r;

            min_i = 1*min(imag(I11(:)));
            I11_i = imag(I11) - min_i;
            scl_i = max(I11_i(:));
            I11_i = I11_i /scl_i;


            [~,I3n_r] = VBM3D((I11_r), sigma_bm3d(kp));% 
            [~,I3n_i] = VBM3D((I11_i), sigma_bm3d(kp));%   


            I3 = double(I3n_r * scl_r + min_r + 1i * (I3n_i * scl_i +  min_i));

    
            tmp_ki = param.E*I3;
            if 1
                I11 = I3 + param.E'*(kspace_raw - tmp_ki*1);
%             
            else
                tmp_tv = param.W*I3;
                I11 = I3 + lambda*tmp_tv + param.E'*(kspace_raw - tmp_ki*1);
            end

        end
      
end

toc

recon_nufft = flipud(recon_nufft);
recon_vBM = flipud(I11);

% display 4 frames
recon_nufft2=recon_nufft(:,:,1);recon_nufft2=cat(2,recon_nufft2,recon_nufft(:,:,7));recon_nufft2=cat(2,recon_nufft2,recon_nufft(:,:,13));recon_nufft2=cat(2,recon_nufft2,recon_nufft(:,:,23));
recon_bm2=recon_vBM(:,:,1);recon_bm2=cat(2,recon_bm2,recon_vBM(:,:,7));recon_bm2=cat(2,recon_bm2,recon_vBM(:,:,13));recon_bm2=cat(2,recon_bm2,recon_vBM(:,:,23));
figure;
subplot(2,1,1),imshow(abs(recon_nufft2),[]);title('Zero-filled FFT')
subplot(2,1,2),imshow(abs(recon_bm2),[]);title('vBM3d')
