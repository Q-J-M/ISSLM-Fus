function continuous_subset = is_consecutive(arr)  
    % 假设arr是一个一维数组  
    % 检查数组是否为空  
    if isempty(arr)  
        error('数组不能为空');  
    end  
      
    % 对数组进行排序，以便更容易检查连续性  
    sorted_arr = sort(arr);  
      
    % 检查数组是否连续  
    consecutive = true;  
    for i = 1:length(sorted_arr)-1  
        if sorted_arr(i+1) ~= sorted_arr(i) + 1  
            consecutive = false;  
            break;  
        end  
    end  
      
    % 如果数组连续，返回从0到末尾的子数组  
    if consecutive  
        continuous_subset = arr;  
    else  
        % 如果不连续，找到断开的位置  
        discontinuity_idx = find(sorted_arr(2:end) ~= sorted_arr(1:end-1) + 1, 1) + 1;  
        % 返回从0到断开位置的子数组  
        continuous_subset = arr(1):sorted_arr(discontinuity_idx-1);  
    end  
end

