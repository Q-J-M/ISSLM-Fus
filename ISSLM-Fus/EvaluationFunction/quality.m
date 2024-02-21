function    [ergas,sam,uiqi,ssim]    =    quality(Z,S,sf)
ratio_ergas = 1.0/sf;

x = Z;
y = S;
% Size, bands, samples 
sz_x = size(x);
n_bands = sz_x(3);
n_samples = sz_x(1)*sz_x(2);

% RMSE
aux = sum(sum((x - y).^2, 1), 2)/n_samples;
rmse_per_band = sqrt(aux);
% rmse_total = sqrt(sum(aux, 3)/n_bands);

% ERGAS
mean_y = sum(sum(y, 1), 2)/n_samples;
ergas = 100*ratio_ergas*sqrt(sum((rmse_per_band ./ mean_y).^2)/n_bands);

% SAM
sam= SpectAngMapper( S, Z );
sam=sam*180/pi;
% num = sum(x .* y, 3);
% den = sqrt(sum(x.^2, 3) .* sum(y.^2, 3));
% sam = sum(sum(acosd(num ./ den)))/(n_samples);

% UIQI - calls the method described in "A Universal Image Quality Index"
% by Zhou Wang and Alan C. Bovik
q_band = zeros(1, n_bands);
for idx1=1:n_bands
    q_band(idx1)=img_qi(S(:,:,idx1), Z(:,:,idx1), 32);
end
uiqi = mean(q_band);
ssim=cal_ssim(S, Z,0,0);
end