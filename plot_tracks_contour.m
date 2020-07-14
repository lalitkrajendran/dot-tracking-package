function plot_tracks_contour(X, Y, U, V, grid_spacing, cmin, cmax)
% Function to display tracked values as a contour plot
%
% INPUTS:
% X, Y: track co-ordinates (pix.)
% U, V: tracked displacements along x and y (pix.)
% grid_spacing: desired spacing between grid points (pix.)
% cmin, cmax: colorbar extents
%
% OUTPUTS:
% None
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % set contour limits
    if ~exist('cmin')
        cmin = 0;
    end

    if ~exist('cmax')
        cmax = 1;
    end

    % calculate contour levels
    contour_levels = linspace(cmin, cmax, 100);

    % interpolate tracks onto grid
    [X_grid, Y_grid, U_grid, V_grid] = interpolate_tracks(X, Y, U, V, grid_spacing);
    
    % display contours
    figure
    contourf(X_grid, Y_grid, sqrt(U_grid.^2 + V_grid.^2), contour_levels, 'edgecolor', 'none')
    colormap(flipud(gray));
    h2 = colorbar;
    caxis([cmin, cmax]);
    set_axes(gca);
    annotate_image(gcf, gca);
    title(h2, '(pix.)')
    title('Displacement');
    
    set(gcf, 'resize', 'off');
    set(gcf, 'Position', [360   584   441   316])
    set(gcf, 'resize', 'off');
    drawnow();
end