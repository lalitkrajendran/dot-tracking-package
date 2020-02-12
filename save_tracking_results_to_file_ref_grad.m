function save_tracking_results_to_file_ref_grad(id_ref, id_grad, id_save_directory, size_ref, size_grad, size_save_directory, tracks, Cp, track_save_directory, save_filename)
    % save identification results
    save(fullfile(id_save_directory, save_filename), 'id_ref', 'id_grad');

    % save sizing results
    save(fullfile(size_save_directory, save_filename), 'size_ref', 'size_grad');

    % save tracking
    save(fullfile(track_save_directory, save_filename), 'tracks', 'Cp');
end