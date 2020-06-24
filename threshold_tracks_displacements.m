function tracks2 = threshold_tracks_displacements(tracks, displacement_threshold)

    % add two columns at the end of the tracks array to hold the
    % validated displacements 
    tracks2 = padarray(tracks, [0, 2], NaN, 'post');
    % copy u displacements
    tracks2(:, 16) = tracks(:, 14);
    % copy v displacements
    tracks2(:, 17) = tracks(:, 15);

    % add another column to hold the validation flag (0 for
    % replaced, 1 for not replaced)
    % tracks2 = tracks;
    tracks2 = padarray(tracks2, [0, 1], 0, 'post');
    % find indices that are outside the specified thershold
    indices = find(abs(tracks2(:, 16)) < displacement_threshold  ...
        & abs(tracks2(:, 17) < displacement_threshold));
        
    tracks2(indices, 16) = NaN;
    tracks2(indices, 17) = NaN;
    
    tracks2(indices, 18) = tracks2(indices, 18) + 1;

end