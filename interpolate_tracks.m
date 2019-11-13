function [X_grid, Y_grid, U_grid, V_grid] = interpolate_tracks(X, Y, U, V, grid_spacing)

    %% calculate min and max x,y co-ordinates of the point cloud
    xmin = ceil(min(X));
    xmax = floor(max(X));

    ymin = ceil(min(Y));
    ymax = floor(max(Y));
    
    %% generate 2D grid
    x_grid = xmin:grid_spacing:xmax;
    y_grid = ymin:grid_spacing:ymax;
    
    [X_grid, Y_grid] = meshgrid(x_grid, y_grid);
    
    %% create interpolants
    F_U = scatteredInterpolant(X, Y, U, 'natural');
    F_V = scatteredInterpolant(X, Y, V, 'natural');
    F_U.ExtrapolationMethod = 'none';
    F_V.ExtrapolationMethod = 'none';
    
    %% interpolate displacements onto the grid
    U_grid = F_U(X_grid, Y_grid);
    V_grid = F_V(X_grid, Y_grid);
    
    %% replace NaNs by linear interpolation
    U_grid = fillmissing(U_grid, 'linear');
    V_grid = fillmissing(V_grid, 'linear');   
end