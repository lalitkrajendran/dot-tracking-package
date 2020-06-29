function dot_props = extract_dot_properties(XYDiameter, dot_index)
% Function to extract dot properties from sizing results
%
% INPUTS:
% XYDiameter: array from sizing results
% dot_index: index of the dot for which the properties are to be extracted
%
% OUTPUTS:
% dot_props: a structure containing the following properties
% X, Y : dot co-ordinates on the image (pix.)
% d_x, d_y: dot diameters along the x and y (pix.)
% I: dot peak intensity (AU)
% R: correlation coefficient - measure of dot tilt
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)

    % dot props
    dot_props = struct;

    % extract x co-ordinates of identified dots
    dot_props.X = XYDiameter(dot_index, 1);
    % extract y co-ordinates of identified dots
    dot_props.Y = XYDiameter(dot_index, 2);
    % extract diameters of identified dots
    dot_props.d_x = XYDiameter(dot_index, 3);
    dot_props.d_y = XYDiameter(dot_index, 4);
    % extract correlation coefficient - measure of dot tilt
    dot_props.R = XYDiameter(dot_index, 5);
    % extract peak intensities of identified dots
    dot_props.I = XYDiameter(dot_index, 6);

end