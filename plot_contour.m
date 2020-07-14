function plot_contour(X, Y, f, grid_spacing, cmin, cmax)
% Function to display tracked values as a contour plot
%
% INPUTS:
% X, Y: track co-ordinates (pix.)
% f: scalar to be plotted as a contour
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

    % create grid
    [X_grid, Y_grid] = create_grid(X, Y, grid_spacing);
    
    % interpolate tracks onto grid
    f_grid = interpolate_to_grid(X, Y, f, X_grid, Y_grid);
    
    % display contours
    figure
    contourf(X_grid, Y_grid, f_grid, contour_levels, 'edgecolor', 'none')
    colormap(flipud(gray));
    h2 = colorbar;
    caxis([cmin, cmax]);
    set_axes(gca);
    annotate_image(gcf, gca);
    
    set(gcf, 'resize', 'off');
    set(gcf, 'Position', [360   584   441   316])
    set(gcf, 'resize', 'off');
    drawnow();
end