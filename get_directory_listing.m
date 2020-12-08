function [files, num_files] = get_directory_listing(filepath, search_string, ignore_string)
% Function to get list of files in a directory matching a search pattern.
%
% INPUTS:
% filepath: directory containing files/folders to be listed
% search_string: string pattern for the file name (e.g. 'a*.png')
% ignore_string: string pattern for files to be ignored (e.g.
% 'staticmask.tif')
%
% OUTPUTS:
% files: array of structures containing file names, paths and size
% num_files: number of files detected in the directory matching the pattern
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)


    % get directory listing
    files = dir(fullfile(filepath, search_string));
    num_files = length(files);

    % remove hidden empty files
    file_counter = 0;
    for file_index = 1:num_files
        if files(file_index).bytes > 0 || files(file_index).isdir
            file_counter = file_counter+1;
            files_temp(file_counter) = copy_struct(files(file_index));
        end        
    end

    if file_counter == 0
        files_temp = [];
    end
    
    % update new set of file names
    files = files_temp;
    num_files = file_counter;

    if nargin > 2
        % remove files containing the ignore string
        file_counter = 0;
        clear files_temp;
        for file_index = 1:num_files
            if ~contains(files(file_index).name, ignore_string)
                file_counter = file_counter+1;
                files_temp(file_counter) = copy_struct(files(file_index));
            end
        end
        if file_counter == 0
            files_temp = [];
        end
        
    % update new set of file names
        files = files_temp;
        num_files = file_counter;
    end
    
    % remove hidden system folders ('.*')
    files = remove_empty_folders(files);
    % update number of files again
    num_files = numel(files);
    fprintf('number of files: %d\n', num_files);

end