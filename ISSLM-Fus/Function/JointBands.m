function [img] = JointBands(bands,band_params)
bandSize = band_params.bandSize;
bandNum = band_params.bandNum;
overlapSize = band_params.overlapSize;
[H, W, S, N] = size(bands);
img_sz = [H, W, bandNum * (bandSize - overlapSize) + overlapSize];
mult = zeros(img_sz);
img = zeros(img_sz);
for i = 1 : bandNum
    m = 1 + (i - 1) * (bandSize - overlapSize);
    mult(:,:,m : m + bandSize - 1) = mult(:,:,m : m + bandSize - 1) + 1;
    img(:,:,m : m + bandSize - 1) = img(:,:,m : m + bandSize - 1) + bands(:,:,:,i);
end
img = img./mult;
end

