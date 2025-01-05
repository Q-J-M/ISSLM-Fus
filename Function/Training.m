function [AR_D,centroids] = Training(LR_HSI,HR_MSI,par)
%% Get some parameter
fprintf('\nI am currently studying spatial information.\n');
fprintf('\nPlease wait a moment...\n');

cls_num = 30; % patch nums
b       = 7; %patch size
ws      = 3; % window size
delta   = par.delta; % threshold of low-frequency information
psf        =    fspecial('gaussian',7,1.5); 
num_blks    =   400*(b^2)*40;
X     =  zeros(0);
Y     =  zeros(0);
M     = par.MN(1);
N     = par.MN(2);
L     = par.L;
l     = size(HR_MSI,3);
%% Get patch
% Get HR_MSI Patch
for i = 1 : l 
    im = HR_MSI(:,:,i);
    [Py,Px] = Get_patches(im, b, psf);
    X       =   [X, Px];% 
    Y       =   [Y, Py];% high-frequency
end
% Get LR_MSI Patch
fprintf(".......");
for i = 1 : L
    im2 = LR_HSI(:,:,i);
    [Py,Px] = Get_patches(im2, b, psf);
    X       =   [X, Px];% 
    Y       =   [Y, Py];
end
% select the smooth blocks
m     =  mean(Y);
d     =  (Y-repmat(m, size(Y,1), 1)).^2;
v     =  sqrt( mean( d ) );
[a, idx]  =  find( v >= delta );
X    =   X(:, idx);
Y    =   Y(:, idx);
%% Clustering.....
[centroids, cls_idx, s_idx, seg, cls_num]   =  SpatialClustering( Y, cls_num );
AR_D     =  zeros(ws*ws-1, cls_num);
%% Auto-regression learning
for  i  =  1 : length(seg)-1
    idx    =   s_idx(seg(i)+1:seg(i+1));
    cls    =   cls_idx(idx(1));
    X1     =   X(:, idx);
    AR    =  Get_AR( X1, ws );
    AR_D(:, cls)    =  AR;
end

