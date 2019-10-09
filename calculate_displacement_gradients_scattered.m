function [X_grid, Y_grid, dU_dx, dU_dy, dV_dx, dV_dy] = calculate_displacement_gradients_scattered(X, Y, U, V, grid_spacing)
% Function to calculate displacement gradients from scattered data
%
% INPUTS:
% X, Y: co-ordinates of the scattered points [pix.]
% U, V: displacements at the scattered points [pix.]
% grid_spacing: desired pixel grid spacing [pix.]
%
% OUTPUTS:
% X_grid, Y_grid: co-ordinates on the grid [pix.]
% dU_dx, dV_dy: displacement gradients [pix./pix.]
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)

    %% interpolate tracks
    [X_grid, Y_grid, U_grid, V_grid] = interpolate_tracks(X, Y, U, V, dot_spacing);
    
    %% calculate displacements gradients
        
    dU_dx = gradient_dudxdvdy_compact_rich_modified(U_grid, grid_spacing, 2);
    dU_dy = gradient_dudxdvdy_compact_rich_modified(U_grid, grid_spacing, 1);

    dV_dx = gradient_dudxdvdy_compact_rich_modified(V_grid, grid_spacing, 2);
    dV_dy = gradient_dudxdvdy_compact_rich_modified(V_grid, grid_spacing, 1);

end