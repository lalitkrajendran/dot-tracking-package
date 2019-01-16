function [x_pos_pix, y_pos_pix] = pos_to_xy_pix_02(x_pos, y_pos, fov_xmin, fov_ymin, pixel_pitch)

    % convert positions to pixels
    x_pos_pix = (x_pos - fov_xmin)/pixel_pitch;
    y_pos_pix = (y_pos - fov_ymin)/pixel_pitch;
        
end