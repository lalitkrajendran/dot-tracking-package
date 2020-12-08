function run_dot_tracking_v5(io, id, sizing, tracking, experimental_parameters)
    % ==========================
    %% check parameter files
    % ==========================    
    [io, id, sizing, tracking] = check_parameter_files(io, id, sizing, tracking);
    
    % ==========================
    %% save parmeters to file
    % ==========================
    % save input/output parameters
    save(fullfile(io.results_save_directory, 'io_params.mat'), 'io');
    
    % save id parameters
    id_save_directory = fullfile(io.results_save_directory, 'id');
    if ~exist(id_save_directory)
        mkdir(id_save_directory);
    end
    save(fullfile(id_save_directory, 'id_params.mat'), 'id');
    
    % save size parameters
    size_save_directory = fullfile(io.results_save_directory, 'size');
    if ~exist(size_save_directory)
        mkdir(size_save_directory);
    end
    save(fullfile(size_save_directory, 'size_params.mat'), 'sizing');
    
    % save track parameters
    track_save_directory = fullfile(io.results_save_directory, 'tracks');
    if ~exist(track_save_directory)
        mkdir(track_save_directory);                
    end
    save(fullfile(track_save_directory, 'track_params.mat'), 'tracking');
    
    % ==========================
    %% extract parameters for id and sizing
    % ==========================    
    % copy settings to data structure
    particleIDprops.v = id.intensity_threshold_current;
    particleIDprops.method = id.segmentation_method; %'dynamic';
    particleIDprops.contrast_ratio = 0;

    % copy settings to data structure
    sizeprops.thresh = id.intensity_threshold_current;
    sizeprops.p_area = id.min_area; %d_p^2;
    sizeprops.sigma = 4;
    sizeprops.errors = double(sizing.default_iwc); % retain IWC estimate if gaussian fit fails
    sizeprops.method = sizing.centroid_subpixel_fit;

    % ==========================
    %% load image mask
    % ==========================
    if io.image_masking
        % load image mask
        image_mask = imread(io.image_mask_filepath);
        % convert mask to 0 and 1
        image_mask(image_mask > 0) = 1;
        % convert mask to double
        image_mask = double(image_mask);
        % flip image mask
        image_mask = flipud(image_mask);
    else
        % create a fake mask that is true everywhere
        image_mask = ones(experimental_parameters.camera_design.y_pixel_number, experimental_parameters.camera_design.x_pixel_number);
    end
    
    % ==========================
    %% calculate dot positions from prior information
    % ==========================    
    % load calibration data
    if strcmp(id.camera_model, 'soloff')
        calibration_results = load(fullfile(id.camera_model_directory, 'calibration_data.mat'));
        experimental_parameters.bos_pattern.X_Min = min(calibration_results.calibration_data.x_world_full{1}) * 1e3;
        experimental_parameters.bos_pattern.Y_Min = min(calibration_results.calibration_data.y_world_full{1}) * 1e3;
    else
        experimental_parameters.bos_pattern.X_Min = 0; %-1.5e3; %0;
        experimental_parameters.bos_pattern.Y_Min = 0; %-3e3; %0;
    end
    
    % calculate expected field of view based on the magnification and the
    % size of the camera sensor
    field_of_view = experimental_parameters.camera_design.x_pixel_number * experimental_parameters.camera_design.pixel_pitch / experimental_parameters.lens_design.magnification;
    experimental_parameters.bos_pattern.dot_number = round(field_of_view/experimental_parameters.bos_pattern.dot_spacing) * 1;
    
    % create array of x,y co-ordinates to describe the dot pattern
    x_array = experimental_parameters.bos_pattern.X_Min + experimental_parameters.bos_pattern.dot_spacing * (0:experimental_parameters.bos_pattern.dot_number);
    y_array = experimental_parameters.bos_pattern.Y_Min + experimental_parameters.bos_pattern.dot_spacing * (0:experimental_parameters.bos_pattern.dot_number);    
    [positions.x, positions.y] = meshgrid(x_array, y_array);
    
    % load camera mapping function coefficients
    if strcmp(id.camera_model, 'soloff')
        mapping_coefficients = load(fullfile(id.camera_model_directory, ['camera_model_type=' num2str(id.order_z) '.mat']));
    else
        mapping_coefficients = [];
    end
    
    % ==========================
    % calculate reference dot locations
    % ==========================
    if io.display_intermediate_progress
        fprintf('calculating reference dot locations\n');
    end
    [pos_ref_dots.x, pos_ref_dots.y] = calculate_reference_dot_locations_new(positions, experimental_parameters, id.camera_model, mapping_coefficients, id.order_z, id.starting_index_x, id.starting_index_y);
    % invert y location of dots due to flipping the image               
    pos_ref_dots.y = experimental_parameters.camera_design.y_pixel_number - pos_ref_dots.y;
        
    % remove points outside FOV
    indices = pos_ref_dots.x < 1 | pos_ref_dots.x > experimental_parameters.camera_design.x_pixel_number-1 | pos_ref_dots.y < 1 | pos_ref_dots.y > experimental_parameters.camera_design.y_pixel_number-1;
    pos_ref_dots.x(indices) = [];
    pos_ref_dots.y(indices) = [];
    num_dots_ref = numel(pos_ref_dots.x);

    % if predicted locations lie in the masked region then ignore
    if io.image_masking
        for dot_index = 1:num_dots_ref
            if image_mask(round(pos_ref_dots.y(dot_index)), round(pos_ref_dots.x(dot_index))) == 0
                pos_ref_dots.x(dot_index) = NaN;
                pos_ref_dots.y(dot_index) = NaN;
            end
        end
        nan_indices = isnan(pos_ref_dots.x) | isnan(pos_ref_dots.y);
        pos_ref_dots.x(nan_indices) = [];
        pos_ref_dots.y(nan_indices) = [];
    end

    % save position of reference dots to structure
    id.x_ref = pos_ref_dots.x;
    id.y_ref = pos_ref_dots.y;

    % ==========================
    %% get directory listings for images, vectors and diameters
    % ==========================    
    % get directory listing for images
    [image_filenames,num_image_files] = get_directory_listing(io.image_directory, 'im*.tif' );
    % calculate number of image pairs to read
    num_image_pairs = num_image_files/2;
    if num_image_pairs > (io.imfend - io.imfstart + io.imfstep)/2
        num_image_pairs = (io.imfend - io.imfstart + io.imfstep)/2;
    end
    % get directory listing for correlation vectors results
    if strcmp(tracking.initialization_method, 'correlation') || strcmp(id.diameter_estimation_method, 'correlation')
        % get directory listing for vectors
        [vector_filenames,~] = get_directory_listing(io.correlation_results.directory, [io.correlation_results.basename '*.mat'], 'corrplane');
    else
        vector_filenames = [];
    end
    
    % ==========================
    %% reference image
    % ==========================    
    if io.display_intermediate_progress
        fprintf('Reference Image\n');
    end

    % image index (the reference image will be the second one in the
    % sequence)
    if strcmp(io.image_type, 'experimental')
        im_index = io.imfstart + 1;
    else
        im_index = io.imfstart;
    end
    
    % load image
    im_ref = imread(fullfile(image_filenames(im_index).folder, image_filenames(im_index).name));
    
    % prepare image for processing
    im_ref = pre_process_image(im_ref, id.minimum_subtraction, id.minimum_intensity_level, image_mask);
        
    % ---------------------------
    %% Identification and Sizing
    % ---------------------------
    if strcmp(id.identification_method, 'standard')
        % ---------------------------------------------
        % standard dot identification and sizing
        % ---------------------------------------------

        % Dot identification
        if io.display_intermediate_progress
            fprintf('running dot identification\n');
        end

        [id_ref.p_matrix, id_ref.peaks, id_ref.num_p] = particle_ID(im_ref,particleIDprops);

        % Dot sizing
        if io.display_intermediate_progress
            fprintf('running dot sizing\n'); 
        end

        % perform dot sizing for frame 1 (only for the first
        % image pair as all reference images are identical)
        [size_ref.XYDiameter, size_ref.mapsizeinfo, size_ref.locxy] = particle_sizing(im_ref, id_ref.p_matrix, id_ref.num_p, sizeprops);

    elseif contains(id.identification_method, 'apriori')
        % ------------------------------------------------------------
        % dot identification and sizing using prior information about
        % the dot locations
        % ------------------------------------------------------------
        
        % create dummy structure for identification results
        id_ref = struct;
        id_ref.x_ref = id.x_ref;
        id_ref.y_ref = id.y_ref;
        
        % predicted dot diameters for the reference image
        d_p = id.dot_diameter*ones(size(pos_ref_dots.x));

        % identify and size dots using their known locations on the
        % target 
        [size_ref.XYDiameter, size_ref.peaks, size_ref.mapsizeinfo, size_ref.locxy, size_ref.mapint] = combined_ID_size_apriori_10(im_ref, pos_ref_dots.x, pos_ref_dots.y, d_p+2, sizing.centroid_subpixel_fit, sizing.default_iwc, id.min_area, id.W_area, id.W_intensity, id.W_distance);
    end
    
    % remove NaN entries
    size_ref = remove_nan_entries_from_sizing(size_ref);
    
    % extract dot co-ordinates
    x_ref = size_ref.XYDiameter(:,1);
    y_ref = size_ref.XYDiameter(:,2);
    z_ref = zeros(size(x_ref));
    
    % ==========================
    %% load and process image pairs
    % ==========================
    
    % number of image pairs to read
    num_image_pairs = (io.imfend - io.imfstart)/io.imfstep + 1;

    for image_pair_index = 1:num_image_pairs        
        % display progress to user
        fprintf('Image Pair Number: %d\n', image_pair_index);
        
        % ==========================
        %% load gradient image
        % ==========================    
        % grad image index
        if strcmp(io.image_type, 'experimental')
            im_index = io.imfstart + 2 * (image_pair_index - 1);
        else
            im_index = io.imfstart + 2 * (image_pair_index - 1) + 1;
        end
        
        % load image
        im_grad = imread(fullfile(image_filenames(im_index).folder, image_filenames(im_index).name));

        % prepare image for processing
        im_grad = pre_process_image(im_grad, id.minimum_subtraction, id.minimum_intensity_level, image_mask);

        % ==========================
        %% load results from correlation processing for the current image pair
        % ==========================
        if strcmp(tracking.initialization_method, 'correlation') || strcmp(id.diameter_estimation_method, 'correlation')
            % load vectors from cross-correlation
            vector_index = (image_pair_index-1)/(io.correlation_frame_step/2) + 1;
            correlation_results_current = load(fullfile(vector_filenames(vector_index).folder, vector_filenames(vector_index).name));
        end

        % ==========================
        %% Identification and Sizing
        % ==========================
        if strcmp(id.identification_method, 'standard')
            % ---------------------------------------------
            % standard dot identification and sizing
            % ---------------------------------------------
            
            %% Dot identification
            
            % perform dot identification for frame 2
            [id_grad.p_matrix, id_grad.peaks, id_grad.num_p] = particle_ID(im_grad, particleIDprops);
                
            %% Dot sizing
            if io.display_intermediate_progress
                fprintf('running dot sizing\n');
            end

            % perform dot sizing for frame 2
            [size_grad.XYDiameter, size_grad.mapsizeinfo, size_grad.locxy]=particle_sizing(im_grad, id_grad.p_matrix, id_grad.num_p, sizeprops);

        elseif contains(id.identification_method, 'apriori')
            % ------------------------------------------------------------
            % dot identification and sizing using prior information about
            % the dot locations
            % ------------------------------------------------------------

            % create dummy structure for identification results
            id_grad = struct;
            
            % estimate effective diameter from correlation
            if strcmp(id.diameter_estimation_method, 'correlation')
                d_p = estimate_effective_dot_diameter(correlation_results_current.X, correlation_results_current.Y, correlation_results_current.Di(:,:,end), pos_ref_dots.x, pos_ref_dots.y, id.dot_diameter);
            else
                d_p = id.dot_diameter*ones(size(pos_ref_dots.x));
            end
            
            % estimate probable locations of dots on frame 2 based on their
            % position in frame 1 and the displacement field from
            % correlation or light ray displacements
            if strcmp(id.initialization_method, 'correlation')
                [x_grad_est, y_grad_est] = calculate_predicted_dot_positions_02(x_ref, y_ref, id.initialization_method, correlation_results_current.X, correlation_results_current.Y, correlation_results_current.U, correlation_results_current.V);
            else
                x_grad_est = x_ref;
                y_grad_est = y_ref;
            end
                       
            % identify and size dots using their known locations on the target
            [size_grad.XYDiameter, size_grad.peaks, size_grad.mapsizeinfo, size_grad.locxy, size_grad.mapint] = combined_ID_size_apriori_10(im_grad, x_grad_est, y_grad_est, d_p+2, sizing.centroid_subpixel_fit, sizing.default_iwc, id.min_area, id.W_area, id.W_intensity, id.W_distance);
        end
        
        % remove NaN entries
        size_grad = remove_nan_entries_from_sizing(size_grad);
        
        % extract co-ordinates
        x_grad = size_grad.XYDiameter(:,1);
        y_grad = size_grad.XYDiameter(:,2);
        z_grad = zeros(size(x_grad));

        % ==========================
        %% Tracking
        % ==========================    
        if io.display_intermediate_progress
            fprintf('running dot tracking\n');
        end

        % if a guess for the dot positions in the second frame is required
        % and does not already exist, then calculae it
        if ~strcmp(tracking.initialization_method, 'none')
            if ~exist('X2_est', 'var') && ~exist('Y2_est', 'var')            
                if strcmp(image_type, 'original')
                    [x_grad_est, y_grad_est] = calculate_predicted_dot_positions(x_ref, y_ref, tracking.initialization_method, correlation_results_current.X, correlation_results_current.Y, correlation_results_current.U, correlation_results_current.V);
                elseif strcmp(image_type, 'deform')
                    x_grad_est = x_ref;
                    y_grad_est = y_ref;
                end
            end
        % if the guess is not required then set it to be the positions of 
        % the dots in the first frame
        else
            x_grad_est = x_ref;
            y_grad_est = y_ref;
        end
        
        % expected z co-ordinates of dots in the second frame (zero)
        z_grad_est = zeros(size(x_ref));

        % extract dot diameters from frame 1
        d_ref = sqrt(size_ref.XYDiameter(:,3).^2 + size_ref.XYDiameter(:,4).^2);
        % extract dot diameters from frame 2
        d_grad = sqrt(size_grad.XYDiameter(:,3).^2 + size_grad.XYDiameter(:,4).^2);
        
        % extract dot intensities from frame 1
        I_ref = size_ref.XYDiameter(:,6);
        % extract dot intensities from frame 2
        I_grad= size_grad.XYDiameter(:,6);

        % track the particles in the image pair using the 3D weighted
        % nearest neighbor tracking method
        [tracks] = weighted_nearest_neighbor3D(x_ref, x_grad, x_grad_est, y_ref, y_grad, y_grad_est,...
            z_ref, z_grad, z_grad_est, d_ref, d_grad, I_ref, I_grad, [tracking.distance_weight, tracking.size_weight, tracking.intensity_weight], tracking.search_radius);

        % ==========================
        %% correct sub-pixel estimate of displacement using a correlation based estimate
        % ==========================    
        % calculate number of tracks
        num_tracks = size(tracks, 1);
        % data structure to hold the correlation plane
        Cp = cell(num_tracks, 1);
        if tracking.perform_correlation_correction
            if io.display_intermediate_progress
                fprintf('Performing Correlation Correction\n');
            end
            % get final sub-pixel displacement estimate by cross-correlation intensity
            % maps of the dots
            U = zeros(num_tracks,1);
            V = zeros(num_tracks,1);
            for track_index = 1:num_tracks
                if rem(track_index, 1000) == 0
                    if io.display_intermediate_progress
                        fprintf('Track: %d of %d\n', track_index, num_tracks);
                    end
                end
                % extract current track
                track_current = tracks(track_index, :);
                [U(track_index), V(track_index), Cp{track_index}] = cross_correlate_dots_07(im_ref, im_grad, size_ref, size_grad, track_current, tracking.correlation_correction.subpixel_fit, tracking.correlation_correction.zero_mean, tracking.correlation_correction.min_sub);            
            end
            % append results to track
            tracks = [tracks, U, V];
        end
        
        % ==========================
        %% flip tracking results to account for reordering reference and gradient images
        % ==========================    
        if strcmp(io.image_type, 'experimental')
            if io.display_intermediate_progress
                fprintf('flipping tracks for experimental data\n');
            end
            % swap x positions
            tracks(:, [1, 2]) = tracks(:, [2, 1]);
            % swap y positions
            tracks(:, [3, 4]) = tracks(:, [4, 3]);
            % swap z positions
            tracks(:, [5, 6]) = tracks(:, [6, 5]);
            % swap diameters
            tracks(:, [7, 8]) = tracks(:, [8, 7]);
            % swap intensities
            tracks(:, [9, 10]) = tracks(:, [10, 9]);
            % swap particle ids
            tracks(:, [11, 12]) = tracks(:, [12, 11]);
            
            if tracking.perform_correlation_correction
                % change sign on displacements from correlation correction
                tracks(:, 14) = -tracks(:, 14);
                tracks(:, 15) = -tracks(:, 15);
            end            
        end
        
        % ==========================
        %% displacement vector validation
        % ==========================            
        % peform validation if specified
        if tracking.perform_validation
            if io.display_intermediate_progress
                fprintf('performing validation\n');
            end
            
            % ==========================
            %% create array to hold results
            % ==========================
            % extract number of columns
            num_cols = size(tracks, 2);
            % add columns at the end of the tracks array to hold the
            % validated displacements 
            tracks = padarray(tracks, [0, 3], NaN, 'post');
            % copy u displacements
            tracks(:, num_cols+1) = tracks(:, num_cols-1);
            % copy v displacements
            tracks(:, num_cols+2) = tracks(:, num_cols);
            % set column to hold the validation flag (0 for % replaced, 1 for not replaced)
            tracks(:, num_cols+3) = 0;
            
            % ==========================
            %% displacement thresholding
            % ==========================
            if io.display_intermediate_progress
                fprintf('performing displacement thresholding\n');
            end

            % add another column to hold the validation flag (0 for
            % replaced, 1 for not replaced)
            tracks = padarray(tracks, [0, 1], 0, 'post');
            if tracking.validation.perform_displacement_thresholding
                if io.display_intermediate_progress
                    fprintf('performing displacement thresholding\n');
                end
                % find indices that are outside the specified thershold
                indices = find(abs(tracks(:, num_cols+1)) > tracking.validation.max_displacement_threshold ...
                    | abs(tracks(:, num_cols+2)) > tracking.validation.max_displacement_threshold);
                
                tracks(indices, num_cols+1) = NaN;
                tracks(indices, num_cols+2) = NaN;
                
                tracks(indices, num_cols+3) = tracks(indices, num_cols+3) + 1;
            end
            
            % ==========================
            %% uod
            % ==========================
            if tracking.validation.perform_uod
                if io.display_intermediate_progress
                    fprintf('performing uod\n');
                end    
                % perform uod
                [u_val, v_val, val, ~, ~] = uod_ptv(tracks(:, 1), tracks(:, 3), tracks(:, num_cols-1), tracks(:, num_cols), ...
                    tracking.validation.replace_vectors, tracking.validation.uod_residual_threshold, tracking.validation.uod_epsilon);
                % update tracks
                tracks(:, num_cols+1) = u_val;
                tracks(:, num_cols+2) = v_val;
                
                if tracking.validation.replace_vectors
                    % add replacement flag
                    tracks(:, num_cols+3) = tracks(:, num_cols+3) + val;
                end
            end            
        end
        
        % ==========================
        %% save results to file
        % ==========================
        if io.display_intermediate_progress
            fprintf('saving results to file\n');
        end
        % name to save the file
        save_filename = ['file_' num2str(im_index, '%04d') '.mat'];
        % save results to file
        save_tracking_results_to_file_ref_grad(id_ref, id_grad, id_save_directory, ...
            size_ref, size_grad, size_save_directory, ...
            tracks, Cp, track_save_directory, save_filename);
        
    end
end
