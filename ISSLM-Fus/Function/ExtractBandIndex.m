function [bandsIndex] = ExtractBandIndex(RM, L)
for i = 1 : L
    band = RM(:, i);
    bandsIndex{i} = find(band > 0.95);
end
end

