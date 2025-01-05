function [bands] = ExtractBand(img,band_params)
bandSize = band_params.bandSize;
bandNum = band_params.bandNum;
overlapSize = band_params.overlapSize;
sz = size(img);
bands = zeros(sz(1),sz(2),bandSize,bandNum);
for i = 1 : bandNum
    m = 1 + (i - 1) * (bandSize - overlapSize);
    bands(:,:,:,i) = img(:,:,m : m + bandSize -1); 
end
end

