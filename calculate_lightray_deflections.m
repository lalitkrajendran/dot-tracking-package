function [ray_displacement, ray_deflection] = calculate_lightray_deflections(pos_1, pos_2, dir_1, dir_2)
% input - filepath containing displacements
% output - deflection in radians

        %% calculate reference displacements from light ray positions
        
        % extract x coordinates (microns)
        pos_1_x = pos_1(1:3:end);
        pos_2_x = pos_2(1:3:end);
        
        % extract y coordinates (microns)
        pos_1_y = pos_1(2:3:end);
        pos_2_y = pos_2(2:3:end);               
        
        % extract z coordinates (microns)
        pos_1_z = pos_1(3:3:end);
        pos_2_z = pos_2(3:3:end);               
        
        ray_displacement = struct;
        
        ray_displacement.x =  (pos_2_x - pos_1_x);
        ray_displacement.y =  (pos_2_y - pos_1_y);

        %% calculate deflections from light ray directions
        
        % extract x direction cosine
        dir_1_x = dir_1(1:3:end);
        dir_2_x = dir_2(1:3:end);
        
        % extract y direction cosine
        dir_1_y = dir_1(2:3:end);
        dir_2_y = dir_2(2:3:end);               
        
        % extract z direction cosine
        dir_1_z = dir_1(3:3:end);
        dir_2_z = dir_2(3:3:end);               
        
        ray_deflection = struct;
        
        ray_deflection.x = abs(acos(dir_2_x) - acos(dir_1_x));
        ray_deflection.y = abs(acos(dir_2_y) - acos(dir_1_y));
        ray_deflection.z = abs(acos(dir_2_z) - acos(dir_1_z));

end