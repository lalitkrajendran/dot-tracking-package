function [X, Y, Z, d, I] = extract_dot_properties(XYDiameter)
% Function to extract dot properties from sizing results
%
% INPUTS:
% XYDiameter: array from sizing results
%
% OUTPUTS:
% X, Y, Z : dot co-ordinates on the image (pix.)
% d: dot diameter (pix.)
% I: dot peak intensity (AU)
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)

    % extract x co-ordinates of identified dots
    X = XYDiameter(:,1);
    % extract y co-ordinates of identified dots
    Y = XYDiameter(:,2);
    % extract z co-ordinates of identified dots
    Z = zeros(size(X));
    % extract diameters of identified dots
    d = XYDiameter(:,3);
    % extract peak intensities of identified dots
    I = XYDiameter(:,4);

end