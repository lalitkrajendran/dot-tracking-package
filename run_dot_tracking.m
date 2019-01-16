function run_dot_tracking(jobfile)
        %% load images


        %% Identification and Sizing

        if strcmp(method, 'iterative')
            [X1, Y1] = identify_dots_iterative(im1, Data, parameters, positions, 'ref', dot_removal_method, threshold_relaxation_parameter, num_dots_ref, d_p/2);
            [X2, Y2] = identify_dots_iterative(im2, Data, parameters, positions, 'grad', dot_removal_method, threshold_relaxation_parameter, num_dots_ref, d_p/2);

            %% remove ghost points
    %             [X1, Y1] = remove_ghost_points(X1, Y1, pos_ref_dots.x, pos_ref_dots.y, tolerance);
            [X1, Y1] = remove_ghost_points(X1, Y1, pos_ref_rays_1.x, pos_ref_rays_1.y, d_p/2);
    %             X1 = X1_iter;
    %             Y1 = Y1_iter;
        elseif strcmp(method, 'standard')
            %% Dot identification
%             particleIDprops.method = 'dynamic';
            fprintf('running dot identification\n');

            [ID1_all{image_pair_index}.p_matrix,ID1_all{image_pair_index}.peaks,ID1_all{image_pair_index}.num_p]=particle_ID(im1,particleIDprops);
            [ID2_all{image_pair_index}.p_matrix,ID2_all{image_pair_index}.peaks,ID2_all{image_pair_index}.num_p]=particle_ID(im2,particleIDprops);

            %% Dot sizing
            fprintf('running dot sizing\n');

            % pos_prana = cell(1,num_sizing_methods);
            pos_prana = [];

            [SIZE1_all{image_pair_index}.XYDiameter,SIZE1_all{image_pair_index}.mapsizeinfo,SIZE1_all{image_pair_index}.locxy]=particle_sizing(im1,ID1_all{image_pair_index}.p_matrix,...
                                ID1_all{image_pair_index}.num_p,sizeprops);
            [SIZE2_all{image_pair_index}.XYDiameter,SIZE2_all{image_pair_index}.mapsizeinfo,SIZE2_all{image_pair_index}.locxy]=particle_sizing(im2,ID2_all{image_pair_index}.p_matrix,...
                                ID2_all{image_pair_index}.num_p,sizeprops);

            % remove nan estimates
            [nan_r, ~] = find(isnan(SIZE1_all{image_pair_index}.XYDiameter));
            SIZE1_all{image_pair_index}.XYDiameter(nan_r, :) = [];
            SIZE1_all{image_pair_index}.mapsizeinfo(nan_r, :) = [];
            SIZE1_all{image_pair_index}.locxy(nan_r, :) = [];

            [nan_r, ~] = find(isnan(SIZE2_all{image_pair_index}.XYDiameter));
            SIZE2_all{image_pair_index}.XYDiameter(nan_r, :) = [];
            SIZE2_all{image_pair_index}.mapsizeinfo(nan_r, :) = [];
            SIZE2_all{image_pair_index}.locxy(nan_r, :) = [];

            X1 = SIZE1_all{image_pair_index}.XYDiameter(:,1);
            Y1 = SIZE1_all{image_pair_index}.XYDiameter(:,2);
            Z1 = zeros(size(X1));
            
            X2 = SIZE2_all{image_pair_index}.XYDiameter(:,1);
            Y2 = SIZE2_all{image_pair_index}.XYDiameter(:,2);
            Z2 = zeros(size(X2));
        else
            if contains(method, 'apriori')
                %% image 1
                
                fprintf('im1\n');

                % estimate effective diameter from correlation
                if strcmp(diameter_estimation_method, 'correlation')
                    d_p = estimate_effective_dot_diameter(diameter_data.X, diameter_data.Y, diameter_data.Di(:,:,end), pos_ref_dots.x, pos_ref_dots.y, d_p_approx);
                else
                    d_p = d_p_approx*ones(size(pos_ref_dots.x));
                end
                
                % minimum area for a set of pixels to be considered a dot
                min_area = 0.5 * median(d_p)^2;
                
                % identify dots using their known locations on the target
                if contains(method, 'rays')
                    [SIZE1_all{image_pair_index}.XYDiameter,SIZE1_all{image_pair_index}.mapsizeinfo,SIZE1_all{image_pair_index}.locxy, SIZE1_all{image_pair_index}.mapint]=combined_ID_size_apriori_03(im1,pos_ref_rays_1.x, pos_ref_rays_1.y, d_p, subpixel_fit, min_area);
                elseif contains(method, 'dots')
%                     [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy, SIZE1_temp{image_pair_index}.mapint]=combined_ID_size_apriori_03(im1,pos_ref_dots.x, pos_ref_dots.y, d_p+2, subpixel_fit, min_area);
%                     [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy, SIZE1_temp{image_pair_index}.mapint]=combined_ID_size_apriori_04(im1,pos_ref_dots.x, pos_ref_dots.y, d_p+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                    [SIZE1_all{image_pair_index}.XYDiameter,SIZE1_all{image_pair_index}.peaks, SIZE1_all{image_pair_index}.mapsizeinfo,SIZE1_all{image_pair_index}.locxy, SIZE1_all{image_pair_index}.mapint]=combined_ID_size_apriori_07(im1,pos_ref_dots.x, pos_ref_dots.y, d_p+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                end

                X1 = SIZE1_all{image_pair_index}.XYDiameter(:,1);
                Y1 = SIZE1_all{image_pair_index}.XYDiameter(:,2);
                %% image 2

                fprintf('im2\n');

                if strcmp(initialization_method, 'rays')                    
                    [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, pos_ref_rays_1.x, pos_ref_rays_1.y, pos_ref_rays_2.x - pos_ref_rays_1.x, pos_ref_rays_2.y - pos_ref_rays_1.y);                
                elseif strcmp(initialization_method, 'correlation')
                    [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, data.X, data.Y, data.U, data.V);
                else
                    X2_est = X1;
                    Y2_est = Y1;
                end
                
                % estimate effective diameter from correlation
                if strcmp(diameter_estimation_method, 'correlation')
                    d_p = estimate_effective_dot_diameter(diameter_data.X, diameter_data.Y, diameter_data.Di(:,:,end), X2_est, Y2_est, d_p_approx);
                else
                    d_p = d_p_approx*ones(size(pos_ref_dots.x));
                end
                                
%                 [SIZE2_temp{image_pair_index}.XYDiameter,SIZE2_temp{image_pair_index}.mapsizeinfo,SIZE2_temp{image_pair_index}.locxy, SIZE2_temp{image_pair_index}.mapint]=combined_ID_size_apriori_03(im2, X2_est, Y2_est, d_p+2, subpixel_fit, min_area);
%                 [SIZE2_temp{image_pair_index}.XYDiameter,SIZE2_temp{image_pair_index}.mapsizeinfo,SIZE2_temp{image_pair_index}.locxy, SIZE2_temp{image_pair_index}.mapint]=combined_ID_size_apriori_04(im2, X2_est, Y2_est, d_p+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                [SIZE2_all{image_pair_index}.XYDiameter,SIZE2_all{image_pair_index}.peaks,SIZE2_all{image_pair_index}.mapsizeinfo,SIZE2_all{image_pair_index}.locxy, SIZE2_all{image_pair_index}.mapint]=combined_ID_size_apriori_07(im2, X2_est, Y2_est, d_p+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                
                X2 = SIZE2_all{image_pair_index}.XYDiameter(:,1);
                Y2 = SIZE2_all{image_pair_index}.XYDiameter(:,2);
            else
                fprintf('unknown input\n');
                keyboard;
            end
        end
        
        % set z co-ordinates to zero
        Z1 = zeros(size(X1));
        Z2 = zeros(size(X2));              
        %% Tracking
        fprintf('running dot tracking\n');
        if ~exist('X2_est') && ~exist('Y2_est')            
            if strcmp(image_type, 'original')
                [X2_est, Y2_est] = calculate_predicted_dot_positions(X1, Y1, initialization_method, data.X, data.Y, data.U, data.V);
            elseif strcmp(image_type, 'deform')
                X2_est = X1;
                Y2_est = Y1;
            end
        end

%         X2_est = X1 + u_guess;
%         Y2_est = Y1 + v_guess;
        Z2_est = Z1;

        d1 = SIZE1_all{image_pair_index}.XYDiameter(:,3);
        d2 = SIZE2_all{image_pair_index}.XYDiameter(:,3);
        I1 = SIZE1_all{image_pair_index}.XYDiameter(:,4);
        I2 = SIZE2_all{image_pair_index}.XYDiameter(:,4);

        weights = [distance_weight, size_weight, intensity_weight]; %[1, 0, 0];

        s_radius = search_radius;

        % run weighted nearest neighbors
        %track the particles in the image pair using the 3D weighted
        %nearest neighbor tracking method
        [tracks_all{image_pair_index}]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,...
            Z1,Z2,Z2_est,d1,d2,I1,I2,weights,s_radius);

        %% correct sub-pixel estimate of displacement using a correlation based estimate
        if correlation_correction
            fprintf('Performing Correlation Correction\n');
            % get final sub-pixel displacement estimate by cross-correlation intensity
            % maps of the dots
            num_tracks = size(tracks_all{image_pair_index}, 1);
            U = zeros(num_tracks,1);
            V = zeros(num_tracks,1);
            for track_index = 1:num_tracks
                if rem(track_index, 1000) == 0
                    fprintf('Track: %d of %d\n', track_index, num_tracks);
                end
    %             fprintf('track number: %d\n', track_index);
                track_current = tracks_all{image_pair_index}(track_index, :);
%                 [U(track_index), V(track_index)] = cross_correlate_dots(im1, im2, track_current(1), track_current(2), track_current(3), track_current(4), track_current(7), track_current(8), correlation_correction_algorithm, subpixel_fit);            
%                 [U(track_index), V(track_index)] = cross_correlate_dots_02(im1, im2, SIZE1_temp{image_pair_index}, SIZE2_temp{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit);            
%                 [U(track_index), V(track_index)] = cross_correlate_dots_03(SIZE1_temp{image_pair_index}, SIZE2_temp{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit);            
%                 [U(track_index), V(track_index)] = cross_correlate_dots_04(SIZE1_temp{image_pair_index}, SIZE2_temp{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit);            
%                 [U(track_index), V(track_index)] = cross_correlate_dots_05(SIZE1_temp{image_pair_index}, SIZE2_temp{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit, zero_mean, min_sub);            
                [U(track_index), V(track_index)] = cross_correlate_dots_06(im1, im2, SIZE1_all{image_pair_index}, SIZE2_all{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit, zero_mean, min_sub);            
            end
            % append results to track
            tracks_all{image_pair_index} = [tracks_all{image_pair_index}, U, V];
        end    
end
