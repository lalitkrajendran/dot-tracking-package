function dot_index = find_dot(X, Y, xmin, xmax, ymin, ymax)
% Function to find the dot with co-ordinates in the desired interval.
%
% INPUTS:
% X, Y: x, y co-ordinates of the dot
% xmin, xmax: x co-ordinate interval
% ymin, ymax: y co-ordinate interval
%   
% OUTPUTS:
% dot_index: index corresponding to the dot in the co-ordinate arrays
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 09/10/19

    % find dot index
    dot_index = find(X > xmin & X < xmax & Y > ymin & Y < ymax);
end