function [pos_1, pos_2, dir_1, dir_2] = load_lightray_data_02(filepath, num_files_read)
    % this function loads light ray data saved at the end of the ray
    % tracing process.
    
    % if the number of files to be read is not specified then read all the
    % files.
    if nargin < 2
        num_files_read = 0;
    end
    %% load light ray positions
    fprintf('loading light ray positions\n');
    % load positions for case w/o density gradients
    pos_1 = read_files(fullfile(filepath, 'light-ray-positions', 'im1'), num_files_read, 'pos*.bin');
    % load positions for case with density gradients
    pos_2 = read_files(fullfile(filepath, 'light-ray-positions', 'im2'), num_files_read, 'pos*.bin');
    
    %% load light ray directions
    fprintf('loading light ray directions\n');
    % load directions for case w/o density gradients
    dir_1 = read_files(fullfile(filepath, 'light-ray-directions', 'im1'), num_files_read, 'dir*.bin');
    % load positions for case with density gradients
    dir_2 = read_files(fullfile(filepath, 'light-ray-directions', 'im2'), num_files_read, 'dir*.bin');

%     % convert direction cosines to angles
%     dir_1 = acos(dir_1);
%     dir_2 = acos(dir_2);
    
end

function data = read_files(filepath, num_files_read, search_string)
    
    % calculate number of light ray position files available
    [files, num_files] = get_directory_listing(filepath, search_string);
    
    if num_files < num_files_read
        % if the number of files to be read is more than the number of
        % files available, then read only the available files.
        fprintf('not enough files to read. reading only %d files\n', num_files);
        num_files_read = num_files;
    elseif num_files_read == 0
        % if the number of files to be read is not specified then read all the
        % files
        num_files_read = num_files;
    end
    
    data = [];
    for file_index = 1:num_files_read
        
        data_temp = read_data_from_file(fullfile(filepath, files(file_index).name));
        data = [data; data_temp];
    end
end

function data = read_data_from_file(filename)
    % open file
    fid = fopen(filename);
    % load data
    data = fread(fid, 'single');
    % close file
    fclose(fid);    
end