function [x_pos, y_pos] = pos_to_xy_pix_03(pos, parameters, starting_index_x, starting_index_y)
    if nargin < 3
        offset_x = 0.5;
        offset_y = 0.5;
        fprintf('no offset specified. setting x offset to %.2f pix and y offset %.2f pix..\n', offset_x, offset_y);
    end
    
    % extract x-y light ray positions separately
    x_u = pos(1:3:end);
    y_u = pos(2:3:end);


    % This is the coordinate of pixel (1,1) [0][0]
    pixel_1_x = -parameters.camera_design.pixel_pitch * (parameters.camera_design.x_pixel_number - 1) / 2.0;
    pixel_1_y = -parameters.camera_design.pixel_pitch * (parameters.camera_design.y_pixel_number - 1) / 2.0;
    
    % Transformation from distorted co-ordinates in image plane (Xd,Yd) to the final image co-ordinates (Xf,Yf)
    x_f = (x_u - pixel_1_x)/parameters.camera_design.pixel_pitch;
    y_f = (y_u - pixel_1_y)/parameters.camera_design.pixel_pitch;

    % reshape into 1D arrays
    x_pos = reshape(x_f, 1, numel(x_f)) - starting_index_x;
    y_pos = reshape(y_f, 1, numel(y_f)) - starting_index_y;
end