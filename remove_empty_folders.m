function folders_new = remove_empty_folders(folders)
% This function removes directories in a listing that are empty such as '.'
% and '._.'
%
% INPUT:
% folders: struct array containing folder information (returned by dir)
%
% OUTPUT:
% folders_new: updatedstruct array without empty folders
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
    
    % calculate number of folders
    num_folders = numel(folders);
    % create new folders struct array
    folders_new = folders;
    % initialize array containing flagged folders for removal
    flag_remove_folders = zeros(1, num_folders);
    
    % loop through folders and flag for removal
    for folder_index = 1:num_folders
        % check if the folder name starts with '.' and flag
        if strcmp(folders(folder_index).name(1), '.')
            flag_remove_folders(folder_index) = 1;
        end
    end
    
    % remove flagged folders
    folders_new(logical(flag_remove_folders)) = [];
end