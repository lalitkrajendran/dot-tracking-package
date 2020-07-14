function [X_grid, Y_grid] = create_grid(X, Y, grid_spacing)
% Function to create a rectangular grid from a set of scattered data points
%
% INPUTS:
% X, Y: scattered co-ordinates (pix.)
% grid_spacing: desired grid spacing (pix.)
% 
% OUTPUTS:
% X_grid, Y_grid: co-ordinates on the rectangular grid (pix.)
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % calculate min and max x,y co-ordinates of the point cloud
    xmin = ceil(min(X));
    xmax = floor(max(X));

    ymin = ceil(min(Y));
    ymax = floor(max(Y));    
    
    % generate 1D co-ordinate arrays for the grid
    x_grid = xmin:grid_spacing:xmax;
    y_grid = ymin:grid_spacing:ymax;
    
    % generate 2D grid
    [X_grid, Y_grid] = meshgrid(x_grid, y_grid);
end