function set_axes(ax)
% This function sets the axis of a figure to be equal and tight. It also 
% ensures that the y axis is pointing upward.
%
% INPUTS:
% ax: handle to the axes of the figure
%
% OUTPUTS:
% None
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    axis(ax, 'equal');
    axis(ax, 'tight');
    set(gca, 'ydir', 'normal');
        
end