function [X, Y, Z, d, I] = nan_dots_outside_ROI(X, Y, Z, d, I, xmin, xmax, ymin, ymax)
% Function to set dots outside the region of interest to NaN
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
    X(indices) = NaN;
    Y(indices) = NaN;
    Z(indices) = NaN;
    d(indices) = NaN;
    I(indices) = NaN;

end


