function f_grid = interpolate_to_grid(X, Y, f, X_grid, Y_grid)
% Function to interpolate a scalar on to a regular grid
%
% INPUTS:
% X, Y: co-ordinates of scattered data
% f: value of the scalar on the scattered data
% X_grid, Y_grid: co-ordinate grid to interpolate data
% 
% OUTPUTS:
% f_grid: scalar interpolated onto the grid
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % create interpolants
    F = scatteredInterpolant(X, Y, f, 'natural');
    F.ExtrapolationMethod = 'none';
    
    % interpolate scalar onto the grid
    f_grid = F(X_grid, Y_grid);
    
    % replace NaNs by linear interpolation
    f_grid = fillmissing(f_grid, 'linear');
end