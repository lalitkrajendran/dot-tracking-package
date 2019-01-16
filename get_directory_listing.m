function [files, num_files] = get_directory_listing(filepath, search_string)
    
    % get directory listing
    files = dir(fullfile(filepath, search_string));
    num_files = length(files);

    % remove hidden empty files
    file_counter = 0;
%     files_temp = [];
    for file_index = 1:num_files
        if files(file_index).bytes > 0 || files(file_index).isdir
            file_counter = file_counter+1;
            files_temp(file_counter) = copy_struct(files(file_index));
        end        
    end

    if file_counter == 0
        files_temp = [];
    end
    
    files = files_temp;
    num_files = file_counter;
    
    fprintf('number of files: %d\n', num_files);

end