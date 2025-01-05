function [img] = JointBands2(Z_bands,bandsIndex,par,interval)
% MN = par.MN;
% M  = MN(1);
% N  = MN(2);
M  = par.m;
N  = par.n;
L  = par.L;
img_sz = [M, N, L];
mult = zeros(img_sz);
img = zeros(img_sz);
n = 1;
% for i = 1 : interval : L
%     mult(:, :, bandsIndex{i}) = mult(:, :, bandsIndex{i}) + 1;
%     img(:, :, bandsIndex{i}) = img(:, :, bandsIndex{i}) + Z_bands{n};
%     n = n + 1;
% end
for i = 1: interval
    mult(:, :, bandsIndex{i}) = mult(:, :, bandsIndex{i}) + 1;
    img(:, :, bandsIndex{i}) = img(:, :, bandsIndex{i}) + Z_bands{n};
    n = n + 1;
end
img = img./mult;
end

