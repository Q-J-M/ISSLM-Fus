function [bandsIndex] = ExtraBandIndex3(RM, L)
arr = zeros(1, L);
flag = true;
i = 1;
j = 1;
start=3;
step=2;
bandsIndex = {};
arrayMap = containers.Map;
while flag
    band = RM(:, i);
    temp = find(band > 0.8);
    sequence =start-step:start-1;
    isInB = ismember(sequence, temp);
    if  all(isInB)
        location=find(temp==start-step);
    else
        location=find(temp==start);
    end
    temp = is_consecutive(temp(location:end));
%     temp = temp(location:end);
    arrayStr1 = num2str(temp);
    key = mlreportgen.utils.hash(arrayStr1);
    if isKey(arrayMap, key)
        isExist=true;
    else
        arrayMap(key) = j;
        isExist=false;
    end
%     isExist=iscontaining(arrayMap,temp);
    if isExist
        start= temp(end)+1;
    else
        start = temp(end);
        bandsIndex{j} = temp;
        j=j+1;
        arr(:,temp) = arr(:,temp)+1;
    end
    i=start;
    flag = any(arr == 0);
end    
end
