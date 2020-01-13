function [X, Y, Z, d, I] = remove_dots_outside_ROI(X, Y, Z, d, I, xmin, xmax, ymin, ymax)
% Function to remove dots outside the region of interest
%
% INPUTS:
% X, Y, Z: dot co-ordinates on the image (pix.)
% d: dot diameter
% I: dot intensity
% xmin, xmax: left and right boundaries (pix.)
% ymin, ymax: bottom and top boundaries (pix.)
%
% OUTPUTS:
% X, Y, Z: dot co-ordinates on the image (pix.)
% d: dot diameter
% I: dot intensity
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)

    % find dots that lie outside the ROI
    indices = X <= xmin | X >= xmax | Y <= ymin | Y >= ymax;

    % remove dot information for this region
    X(indices) = [];
    Y(indices) = [];
    Z(indices) = [];
    d(indices) = [];
    I(indices) = [];

end


