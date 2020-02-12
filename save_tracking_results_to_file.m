function save_tracking_results_to_file(id1, id2, id_save_directory, size1, size2, size_save_directory, tracks, track_save_directory, save_filename)
    % save identification results
    save(fullfile(id_save_directory, save_filename), 'id1', 'id2');

    % save sizing results
    save(fullfile(size_save_directory, save_filename), 'size1', 'size2');

    % save tracking
    save(fullfile(track_save_directory, save_filename), 'tracks');
end