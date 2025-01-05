function [Z] = ISSLM(LR_HSI,HR_MSI,par,ImgOri,NAR)
%% GET parameter
MN = par.MN;
M  = MN(1);
N  = MN(2);
L  = par.L;
K  = par.K;
K2 = par.K2;
D  = par.D;
R  = par.R;
sf = par.sf;
mu = 0.0001;
lambda1 = par.lambda1;
lambda2 = par.lambda2;
eta2 = par.eta2;
eta1 = par.eta1;
ro = par.ro;
%% Group MSI by non-local similarity
% %Group parameters
patchsize = 8;
overlap = 4;
patch_params.block_sz = [patchsize,patchsize];
patch_params.overlap_sz = [overlap,overlap];
aug1 = (M - patchsize) / (patchsize - overlap) + 1;
aug2 = (N - patchsize) / (patchsize - overlap) + 1;
patch_params.block_num = [aug1,aug2];
% Group patch by kmeans++
HR_MSI3D = hyperConvert3D(HR_MSI, M, N, L);
fkmeans_omaxpt.careful = 1;
preconstruc_blocks = ExtractBlocks(HR_MSI3D,patch_params);
MSI_Patchs = Unfold(preconstruc_blocks,size(preconstruc_blocks),4);
[group] = fkmeans(MSI_Patchs,K2,fkmeans_omaxpt);
%% Group HSI by divide sub-bands
fprintf('\nI am currently studying spectral information.\n');
RM = corrcoef(LR_HSI');
bandsIndex = ExtraBandIndex3(RM,L); 
% interval = 2;
interval=size(bandsIndex,2);
%% initial value
LR_HSI3D = hyperConvert3D(LR_HSI, M / sf, N /sf, L);
HR_load1  =   imresize(LR_HSI3D, sf,'bicubic');
ZE        =   hyperConvert2D(HR_load1);
Z         =   ZE;
W         =   ZE;
A = zeros(K , M * N);
S = zeros(K , M * N);
D0 = R * D;
D02 = D0'*D0;
D2 = D'*D;
Ek = eye(K);
XHT = par.HT(LR_HSI);
D0TY = D0' * HR_MSI;
DTU = D' * ZE;
V1 = zeros( size(Z) );
V2 = zeros( size(A) );
V3 = zeros( size(Z) );
fprintf('\nCommencing model solving...\n');
fprintf('\nIt may take some time, please be patient...\n');
%% iteration
for i = 1:par.Iter2
    A = inv( D02  +   (lambda1 + mu)*D2   + mu*Ek ) * ( D0TY  + lambda1*DTU  + D'*(mu*Z-V1/2) + (mu*S+V2/2) );
    DA    =    D*A;
    Zt    =    hyperConvert3D(DA, M, N, K);
    rmse2(i)=getrmse(double(im2uint8(ImgOri)),double(im2uint8(Zt))); 
    disp(rmse2(i));
    for  j  =  1 : L
        filename = ['G:\paper_code\Lib\',num2str(j),'.mat'];
        file = load(filename);
        Mi = file.ARM;
        MZ = Mi * DA(j,:)';
        MZ = MZ';
        B     =    (XHT(j,:) + mu*(DA(j,:) + V1(j,:) / (2*mu)) + mu * (W(j,:) + V3(j,:) / (2*mu)) + lambda2 * MZ )';% + mu * (W(j,:) + V3(j,:) / (2*mu))  + lambda2 * MZ 
        [z,flag]     =    pcg( @(x)A_x(x, mu, par.fft_B, par.fft_BT, sf, MN, lambda2), B, 1E-3, 350, [], [], Z(j, :)' );
        Z(j, :)      =    z';
    end
    Z(Z<0) = 0;
    Z(Z>1) = 1;
    S = max( soft(A - V2 / (2*mu), eta2/(2*mu)), 0);
%********************************************************************************************
% % fourth method : use TNN
    W = Z -  V3 / (2*mu);
    W = hyperConvert3D(W, M, N, L);
    preconstruc_blocks1 = ExtractBlocks(W,patch_params);
    Z_block1 = zeros(patch_params.block_sz(1),patch_params.block_sz(2),L,patch_params.block_num(1) * patch_params.block_num(2));
    for k = 1 : max(group)
        group_idx = find(group == k);
        Temp = preconstruc_blocks1(:,:,:,group_idx);
        [a,b,c,d] = size(Temp);
        Temp = reshape(Temp,[a*b,c,d]);
        par.m = d;
        par.n = a*b;
        Temp = permute(Temp,[1,3,2]);
        bands = ExtractBand2(Temp, bandsIndex, L, interval);
        nums = max(size(bands));
        for n = 1 : nums
            temp = bands{n};
            temp = permute(temp,[3,2,1]);
            V_LR4 = Log_prox_tnn( temp, eta1/2/mu );
            V_SLR = permute(V_LR4,[2,3,1]);
            Z_bands{n} = V_SLR;
        end
        V_LR1 = JointBands2(Z_bands,bandsIndex, par, interval);
        V_LR = permute(V_LR1,[2,3,1]);
        V_LR = reshape(V_LR,[a,b,c,d]);
        Z_block1(:,:,:,group_idx) = V_LR;
    end
    W = JointBlocks(Z_block1,patch_params);

    W = hyperConvert2D(W);
    V2      =    V2 + mu*( S - A );    
    V1      =    V1 + mu*( DA - Z);
    V3      =    V3 + mu*( W - Z);
    U       =    Z * NAR;
    DTU     =    D'*U;
    mu      =    mu*ro;
end
Z = Zt;
end

