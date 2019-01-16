function d_p = estimate_effective_dot_diameter(X, Y, D, x_p, y_p, d_default)
% This function estimates the effective diameter of a dot in the image
% using the diameter of the correlation plane

% INPUTS:
% X, Y - X and Y location of vectors from correlation [pix.]
% D - Diameter of the correlation peak [pix.]
% x_p, y_p - dot locations in the image (can be approximate) [pix.]
% d_default - default diameter to be assigned in case the interpolation
% fails [pix.]

    % estimate dot diameter from correlation peak by interpolation
    d_p = 1/sqrt(2) * interp2(X, Y, D, x_p, y_p, 'spline');
    
    % set nan elements to be the default diameter specified by the user
    d_p(isnan(d_p)) = d_default;    

end