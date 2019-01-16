function pos_sorted = sort_positions(pos)
    % convert positions to column arrays
    if size(pos.x, 1) == 1
        pos.x = pos.x';
        pos.y = pos.y';
    end
%     pos_array = [pos.x, pos.y];
%     pos_array_sorted = sortrows(pos_array, 'ascend');
%     pos_sorted.x = pos_array_sorted(:,1);
%     pos_sorted.y = pos_array_sorted(:,2);

    % sort the particle positions identified by prana
    [~,idx] = sort(pos.x);
    pos_sorted.x = pos.x(idx);                        
    pos_sorted.y = pos.y(idx);
  
    
end