function [x_c, y_c] = calculate_centroids_from_lightrays_06(x_pos, y_pos, num_lightrays_per_dot) 
% This functions takes in a set of light ray positions on the camera sensor
% and calculates the centroid of the dots by grouping together light rays
% that are within a distance threshold
%
% INPUTS:
% x_pos, y_pos - light ray positions. column arrays. 1x number of light rays
% num_lightrays_per_dot - number of light rays per dot in the image
% generation process. 
%
% OUTPUTS:
% x_c, y_c - co-ordinates of dot centroids
%

     
    % total number of light rays
    num_lightrays_total = length(x_pos);

    % calculate number of dots
    num_dots = num_lightrays_total/num_lightrays_per_dot;
    
    fprintf('number of dots identified from light rays: %d\n', num_dots);
    
    % arrays to hold dot centroids
    x_c = zeros(num_dots,1);
    y_c = zeros(num_dots,1);
   
    
    for dot_index = 1:num_dots
        lightray_start_read_index = (dot_index-1) * num_lightrays_per_dot + 1;
        lightray_stop_read_index = dot_index * num_lightrays_per_dot;
                
        % extract x,y co-ordinates for the current dot
        x_lightrays = x_pos(lightray_start_read_index : lightray_stop_read_index);
        y_lightrays = y_pos(lightray_start_read_index : lightray_stop_read_index);

        % calculate x,y centroids for the current dot without outliers
        x_c(dot_index) = nanmean(x_lightrays);
        y_c(dot_index) = nanmean(y_lightrays);
        
    end


end