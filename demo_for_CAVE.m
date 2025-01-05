clear
clc
addpath(genpath('Function'));
addpath(genpath('EvaluationFunction'));
%% data read
ImgOriFile = load('.\Data\CAVA\sponges.mat');
ImgOri = ImgOriFile.S;
[M,N,L] = size(ImgOri);
MN = [M,N];
sf = 8;% downsampling factor
s0 = 2;
%*********************************************************************************
% some Parameter (parameter may need to regulation add (*))
par.delta = 0.006;%  (*) 
par.lambda = 0.001; % (*)
par.K = 80; % (*) dictionary num
par.K2 = 100;% non-local cluster num
par.eta2 = 0.000006;% sparse coefficients (*)
par.lambda1    =   0.015;% NonLocal-Autoregression item coefficients (*)
par.lambda2    =   0.0008;% Local-Autoregression item coefficients (*)
par.eta1       =   0.01;
par.ro         =   1.1;
par.Iter2       =    20;
par.c2         =    1;
% low-dimension
R = create_F();
%***********************************************************************************
% Blur kernel
psf        =    ones(sf)/(sf^2);
par.psf = psf;
par.fft_B = psf2otf(psf,MN);
par.fft_BT = conj(par.fft_B);
par.H = @(z)H_z(z,par.fft_B,sf,MN,s0);
par.HT = @(y)HT_y(y, par.fft_BT, sf, MN);
par.R = R;
par.MN = MN;
par.L = L;
par.sf = sf;
%% Data simulation
%LR-HSI
ImgOri2D = hyperConvert2D(ImgOri);
LR_HSI = par.H(ImgOri2D);
%HR-MSI
HR_MSI = R * ImgOri2D;
%% Data Training
t0=clock;
LR_HSI3D = hyperConvert3D(LR_HSI,M/sf,N/sf);
HR_MSI3D = hyperConvert3D(HR_MSI, M, N);
[AR_D,centroids] = Training(LR_HSI3D, HR_MSI3D, par);
par.AR_D = AR_D;
par.centroids = centroids;
HR_load1 = imresize(LR_HSI3D, sf,'bicubic');
% Auto-regression Matrix learning
for i = 1 : L
  ARM           =    Compute_AR_Matrix(HR_load1(:,:,i),AR_D,centroids,par);
  filename = ['E:\MyCode2\Lib\',num2str(i),'.mat'];
  save(filename,'ARM');
end
%% NonLocal Auto-regression Matrix learning
NAR = Comp_NLAR_Matrix(HR_MSI,MN);
%% Dictionary Learning
D = Nonnegative_DL(LR_HSI, par);
par.D = D;
%% Use the fuing method
Z = ISSLM(LR_HSI,HR_MSI,par,ImgOri,NAR);
t1=etime(clock,t0);
[psnr4,rmse4, ergas4, sam4, uiqi4,ssim4] = quality_assessment(double(im2uint8(ImgOri)), double(im2uint8(Z)), 0, 1.0/sf);

