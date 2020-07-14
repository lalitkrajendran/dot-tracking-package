function [X_grid, Y_grid, U_grid, V_grid] = interpolate_tracks(X, Y, U, V, grid_spacing)
% Function to interpolate tracked displacements onto a grid
%
% INPUTS:
% X, Y: x, y co-ordinates of the tracks (pix.)
% U, V: x, y displacements of the tracks (pix.)
% grid_spacing: desired grid spacing - usually the dot spacing (pix.)
% 
% OUTPUTS:
% X_grid, Y_grid: co-ordinate grid
% U_grid, V_grid: interpolated displacements on the grid
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % create grid
    [X_grid, Y_grid] = create_grid(X, Y, dot_spacing)

    % interpolate displacements
    U_grid = interpolate_to_grid(X, Y, U, X_grid, Y_grid);
    V_grid = interpolate_to_grid(X, Y, V, X_grid, Y_grid);

end