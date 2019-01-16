function [x, y] = identify_dots_iterative(im, Data, parameters, positions, image_type, dot_removal_method, relaxation_parameter, num_dots_ref, tolerance)
%     num_dots_ref = parameters_ref.bos_pattern.grid_point_number;
%     num_dots_ref = parameters.bos_pattern.x_grid_point_number * parameters.bos_pattern.y_grid_point_number;
    d_diff = parameters.camera_design.diffraction_diameter;
    particleIDprops = Data.ID;
    sizeprops = Data.Size;
   
    %% iterative dot identification
    loop_ctr = 0;
    loop_ctr_max = 100;
    num_dots_identified = 0;
    num_dots_prev = 0;
    num_dots_new = 1;
    im_0 = im;
    size_results_all = [];

    while loop_ctr < loop_ctr_max && num_dots_new > 0 % num_dots_identified < num_dots_ref && num_dots_new > 0
        loop_ctr = loop_ctr + 1;
        % if there are no particles left, then exit
        if sum(im>Data.ID.v) == 0
            fprintf('no dots left in the image. exiting.\n');
            break;
        end

        %% perform identification

        fprintf('Performing Dot Identification ...\n');
        % perform identification
        [ID.p_matrix,ID.peaks,ID.num_p]=particle_ID(im,particleIDprops);

        %% perform dot sizing

        fprintf('Performing Dot Sizing ...\n');
        % run sizing
        [SIZE.XYDiameter,SIZE.mapsizeinfo,SIZE.locxy]=particle_sizing(im, ID.p_matrix, ID.num_p, sizeprops);
        % remove dots that are smaller than the known size
        indices = isnan(SIZE.XYDiameter(:,3)) | SIZE.XYDiameter(:,3) < Data.Size.size_threshold_fraction * d_diff;
        SIZE.XYDiameter(indices,:) = [];

        % aggregate positions
        if loop_ctr == 1
            size_results_all = [size_results_all; SIZE.XYDiameter];
        else
            xy_old = size_results_all(:,1:2);
            xy_new = SIZE.XYDiameter(:,1:2);
            
            [tracks]=weighted_nearest_neighbor3D(xy_new(:,1), xy_old(:,1), xy_new(:,1), ...
                xy_new(:,2), xy_old(:,2), xy_new(:,2), ...
                zeros(size(xy_new(:,1))), zeros(size(xy_old(:,1))), zeros(size(xy_new(:,1))), ...
                  zeros(size(xy_new(:,1))), zeros(size(xy_old(:,1))),...
                zeros(size(xy_new(:,1))), zeros(size(xy_old(:,1))),...
                [1,0,0], tolerance);
            
            % find unmatched particles
            unmatched_indices = find_unmatched_particles(tracks, size(xy_new, 1));
            
            size_results_all = [size_results_all; SIZE.XYDiameter(unmatched_indices,:)];
        end
        
        % number of identified dots
        num_dots_identified = size(size_results_all, 1);
        num_dots_new = num_dots_identified - num_dots_prev;
        num_dots_prev = num_dots_identified;
        fprintf('loop: %d, Dots identified: %d out of %d\n', loop_ctr, num_dots_identified, num_dots_ref);

        %% check if all dots have been identified and remove identified dots

        if dot_removal_method == 1
            % set pixels identified as particles to zero
            particles = ID.p_matrix > 0;
            im_new = im;
            im_new(particles) = 0;

        elseif dot_removal_method == 2
            peak_intensity = max(im(:));
            
            % generate an image formed by the identified dots alone
%             evalc('im_temp = generateParticleImage(size(im,1), size(im,2), SIZE.XYDiameter(:,1), SIZE.XYDiameter(:,2), d_diff*ones(size(SIZE.XYDiameter(:,1))), SIZE.XYDiameter(:,4));');                                    
%             im_temp = generateParticleImage(size(im,1), size(im,2), SIZE.XYDiameter(:,1), SIZE.XYDiameter(:,2), d_diff*ones(size(SIZE.XYDiameter(:,1))), SIZE.XYDiameter(:,4)/0.8);                                  
%             im_temp = generateParticleImage(size(im,1), size(im,2), SIZE.XYDiameter(:,1), SIZE.XYDiameter(:,2), SIZE.XYDiameter(:,3), SIZE.XYDiameter(:,4));                                  
            im_temp = generateParticleImage(size(im,1), size(im,2), size_results_all(:,1), size_results_all(:,2), size_results_all(:,3), size_results_all(:,4));                                  

            % calculate the residual image once the identified dots have been
            % removed
            im_residual_d = double(im_0) - im_temp;
            % remove zero intensity values
            im_residual_d(im_residual_d < 0) = 0;
            % convert to uint16
            im_new = cast(im_residual_d, 'like', im);
        end

        % update the image
        im = im_new;
        
        % reduce threshold so more dots can be identified
        Data.ID.v = Data.ID.v * (1 - relaxation_parameter);
        Data.Size.thresh   = Data.Size.thresh * (1 - relaxation_parameter);

    end
    
    x = size_results_all(:,1);
    y = size_results_all(:,2);
    
end

        % remove repeated entries
%         size_results_all = uniquetol(size_results_all, 0.001/max(abs(size_results_all(:))), 'byrows', true);
%         size_results_all = uniquetol(size_results_all, 0.005, 'byrows', true);
%         size_results_all = unique(size_results_all, 'rows');
