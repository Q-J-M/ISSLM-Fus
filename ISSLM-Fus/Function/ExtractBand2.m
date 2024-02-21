function [bands] = ExtractBand2(W, bandsIndex, L, interval)
n = 1;
for i = 1 : interval : L
    Index = bandsIndex{i};
    bands{n} = W(:,:,Index);
    n = n+1;
%     bands{i} = W(:,Index,:);
end
end

