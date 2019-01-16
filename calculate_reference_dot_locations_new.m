function [x_c, y_c] = calculate_reference_dot_locations_new(positions, parameters, camera_model, mapping_coefficients, orderz, starting_index_x, starting_index_y)

    if ~exist('starting_index_x') || ~exist('starting_index_x')
        starting_index_x = 0;
        starting_index_y = 0;
    end
    if strcmp(camera_model, 'thin-lens')

        % -----------------------------------------------------------------
        % THIN LENS MODEL
        % -----------------------------------------------------------------
        
        % size of a pixel on the camera (um)
        pixel_pitch = parameters.camera_design.pixel_pitch;

        % image field of view (um)
        fov_xmin = parameters.bos_pattern.X_Min; % + pixel_pitch;
        fov_ymin = parameters.bos_pattern.Y_Min; % + pixel_pitch;    

        % This is the coordinate of pixel (1,1) [0][0]
        pixel_1_x = -parameters.camera_design.pixel_pitch * (parameters.camera_design.x_pixel_number - 1) / 2.0;
        pixel_1_y = -parameters.camera_design.pixel_pitch * (parameters.camera_design.y_pixel_number - 1) / 2.0;

        % calculate magnification
        if isfield(parameters.lens_design, 'magnification')
            M = parameters.lens_design.magnification;
        else
            M = parameters.lens_design.focal_length/(parameters.lens_design.object_distance - parameters.lens_design.focal_length);
        end

        % convert positions to x,y arrays in pixels
        [x_c, y_c] = pos_to_xy_pix_02(positions.x, positions.y, fov_xmin, fov_ymin, pixel_pitch);

        % account for magnification
        x_c = M * x_c; % - 1.5;
        y_c = M * y_c; % - 1.5;
        % reshape into 1D arrays
        x_c = reshape(x_c, 1, numel(x_c));
        y_c = reshape(y_c, 1, numel(y_c));

        % invert positions since the camera forms a real image
        x_c = parameters.camera_design.x_pixel_number - x_c;
        y_c = parameters.camera_design.y_pixel_number - y_c;

        x_c = x_c - starting_index_x;
        y_c = y_c - starting_index_y;
    elseif strcmp(camera_model, 'tsai')

        % size of a pixel on the camera (um)
        pixel_pitch = parameters.camera_design.pixel_pitch;

        % image field of view (um)
        fov_xmin = parameters.bos_pattern.X_Min; % + pixel_pitch;
        fov_ymin = parameters.bos_pattern.Y_Min; % + pixel_pitch;    

        % This is the coordinate of pixel (1,1) [0][0]
        pixel_1_x = -parameters.camera_design.pixel_pitch * (parameters.camera_design.x_pixel_number - 1) / 2.0;
        pixel_1_y = -parameters.camera_design.pixel_pitch * (parameters.camera_design.y_pixel_number - 1) / 2.0;

        % -----------------------------------------------------------------
        % TSAI CAMERA MODEL
        % -----------------------------------------------------------------

        % calculate the 3d co-ordinates of the dots in the image co-ordinates
        [X_i,Y_i,Z_i]=rotate_coordinates(positions.x,positions.y,parameters.lens_design.object_distance*ones(size(positions.x)), ...
            parameters.camera_design.x_camera_angle, parameters.camera_design.y_camera_angle, 0, 0, 0, 0);

        image_distance = (1.0/parameters.lens_design.focal_length - 1.0/parameters.lens_design.object_distance)^-1;   
        % Transformation from 3D image co-ordinates to undistorted image plane (Xu,Yu) co-ordinates
        x_u = -image_distance * X_i./Z_i;
        y_u = -image_distance * Y_i./Z_i;

        % Transformation from distorted co-ordinates in image plane (Xd,Yd) to the final image co-ordinates (Xf,Yf)
        x_f = (x_u - pixel_1_x)/parameters.camera_design.pixel_pitch;
        y_f = (y_u - pixel_1_y)/parameters.camera_design.pixel_pitch;

        % reshape into 1D arrays
        x_c = reshape(x_f, 1, numel(x_f)) - starting_index_x;
        y_c = reshape(y_f, 1, numel(y_f)) - starting_index_y;

    elseif strcmp(camera_model, 'soloff')

        % -----------------------------------------------------------------
        % SOLOFF MODEL
        % -----------------------------------------------------------------

%         % calculate the 3d co-ordinates of the dots in the image co-ordinates
%         [x_i,y_i,z_i]=rotate_coordinates(positions.x,positions.y,parameters.lens_design.object_distance*ones(size(positions.x)), ...
%             parameters.camera_design.x_camera_angle, parameters.camera_design.y_camera_angle, 0, 0, 0, 0);
% 
%         % reshape into 1D arrays and convert to mm
%         x_i = reshape(x_i, 1, numel(x_i))/1e3;
%         y_i = reshape(y_i, 1, numel(y_i))/1e3;
%         z_i = zeros(size(x_i));
% 
%         % -----------------------------------------------------------------
%         % load camera mapping function
%         % -----------------------------------------------------------------
%         x1 = x_i;
%         x2 = y_i;
%         x3 = z_i;
        
        x1 = reshape(positions.x, 1, numel(positions.x))/1e3;
        x2 = reshape(positions.y, 1, numel(positions.y))/1e3;
        x3 = zeros(size(x1));
        
        % calculate x,y location of dots on the camera sensor
        if orderz==1                % cubic xy, linear z
            a = mapping_coefficients.aXcam1;
            x_c =a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
            a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
            a(10).*x1.^3 + a(11).*x1.^2.*x2 + a(12).*x1.*x2.^2 +...
            a(13).*x2.^3 + a(14).*x1.^2.*x3 + a(15).*x1.*x2.*x3 +...
            a(16).*x2.^2.*x3;

            a = mapping_coefficients.aYcam1;
            y_c =a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
            a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
            a(10).*x1.^3 + a(11).*x1.^2.*x2 + a(12).*x1.*x2.^2 +...
            a(13).*x2.^3 + a(14).*x1.^2.*x3 + a(15).*x1.*x2.*x3 +...
            a(16).*x2.^2.*x3;

        elseif orderz==2            % cubic xy, quadratic z
            a = mapping_coefficients.aXcam1;
            x_c =a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
                a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
                a(10).*x3.^2 + a(11).*x1.^3 + a(12).*x1.^2.*x2 + a(13).*x1.*x2.^2 +...
                a(14).*x2.^3 + a(15).*x1.^2.*x3 + a(16).*x1.*x2.*x3 +...
                a(17).*x2.^2.*x3 + a(18).*x1.*x3.^2 + a(19).*x2.*x3.^2;    

            a = mapping_coefficients.aYcam1;
            y_c =a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
            a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
            a(10).*x3.^2 + a(11).*x1.^3 + a(12).*x1.^2.*x2 + a(13).*x1.*x2.^2 +...
            a(14).*x2.^3 + a(15).*x1.^2.*x3 + a(16).*x1.*x2.*x3 +...
            a(17).*x2.^2.*x3 + a(18).*x1.*x3.^2 + a(19).*x2.*x3.^2;    
        end

        % flip co-ordinates
%         x_c = parameters.camera_design.x_pixel_number - x_c;
%         y_c = parameters.camera_design.y_pixel_number - y_c;

    else
        fprintf('unknown camera model. exiting.\n');
        return;
    end
end