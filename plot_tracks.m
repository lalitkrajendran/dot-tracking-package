function plot_tracks(tracks, scale_factor, plot_corrected_tracks, plot_validated_tracks)
% Function to plot tracks
%
% INPUTS:
% tracks: result from weighted_nearest_neighbor3D
% scale_factor: factor to multiply displacements. Default = 1.
% plot_corrected_tracks: boolean to specify if the correlation
% corrected tracks (in the 14 and 15 columns are to be plotted).
% plot_validated_tracks: boolean to specify if the validated results are to
% be plotted (in columns 16 and 17)
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 09/10/19

    if nargin < 2
        scale_factor = 1;
        plot_corrected_tracks = false;
        plot_validated_tracks = false;
    elseif nargin < 3
        plot_corrected_tracks = false;
        plot_validated_tracks = false;        
    elseif nargin < 4
        plot_validated_tracks = false;
    end
    
    % extract co-ordinates
    X_track = tracks(:, 1);
    Y_track = tracks(:, 3);
    
    % calculate displacements
    if plot_corrected_tracks
        U_track = tracks(:, 14);
        V_track = tracks(:, 15);
    else
        U_track = (tracks(:, 2) - tracks(:, 1));
        V_track = (tracks(:, 4) - tracks(:, 3));
    end
    
    if size(tracks, 2) >= 16 && plot_validated_tracks
        U_track = tracks(:, 16);
        V_track = tracks(:, 17);
    end
    
    % plot displacements    
    quiver(X_track, Y_track, U_track * scale_factor, V_track * scale_factor, 'autoscale', 'off')
    annotate_image(gcf, gca);
    
    h = gcf;
end