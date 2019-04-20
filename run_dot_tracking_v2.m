function run_dot_tracking_v2(io, id, sizing, tracking, experimental_parameters)
    %% get directory listings for images, vectors and diameters

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
    end

    if io.image_masking
        % load image mask
        image_mask = imread(io.image_mask_filepath);
        % convert mask to 0 and 1
        image_mask(image_mask > 0) = 1;
        % convert mask to double
        image_mask = double(image_mask);
    end
    
    save(fullfile(io.results_save_directory, 'io_params.mat'), 'io');
    %create write directory if it does not exist
    id_save_directory = fullfile(io.results_save_directory, 'id');
    mkdir(id_save_directory);
    save(fullfile(id_save_directory, 'id_params.mat'), 'id');
    
    size_save_directory = fullfile(io.results_save_directory, 'size');
    mkdir(size_save_directory);
    save(fullfile(size_save_directory, 'size_params.mat'), 'sizing');
    
    track_save_directory = fullfile(io.results_save_directory, 'tracks');
    mkdir(track_save_directory);                
    save(fullfile(track_save_directory, 'track_params.mat'), 'tracking');
    
    %% calculate dot positions from prior information
    
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
    
    % calculate reference dot locations
    fprintf('calculating reference dot locations\n');
    [pos_ref_dots.x, pos_ref_dots.y] = calculate_reference_dot_locations_new(positions, experimental_parameters, id.camera_model, mapping_coefficients, id.order_z, id.starting_index_x, id.starting_index_y);
        
    % remove points outside FOV
    indices = pos_ref_dots.x < 1 | pos_ref_dots.x > experimental_parameters.camera_design.x_pixel_number-1 | pos_ref_dots.y < 1 | pos_ref_dots.y > experimental_parameters.camera_design.y_pixel_number-1;
    pos_ref_dots.x(indices) = [];
    pos_ref_dots.y(indices) = [];
    num_dots_ref = numel(pos_ref_dots.x); % - sum(indices);
        
    id.x_ref = pos_ref_dots.x;
    id.y_ref = pos_ref_dots.y;

    %% load and process image pairs

    % initialize data structures to hold results
    ID1_all = cell(1,num_image_pairs);
    ID2_all = cell(1,num_image_pairs);
    SIZE1_all = cell(1,num_image_pairs);
    SIZE2_all = cell(1,num_image_pairs);
    tracks_all = cell(1,num_image_pairs);
    image_pair_index = 0;
    for im_index = io.imfstart:io.imfstep:io.imfend %num_image_pairs        
        image_pair_index = image_pair_index + 1;
        % display progress to user
        fprintf('Image Pair Number: %d\n', image_pair_index);
        
        %% load image pairs and perform masking

        % load images
        im1 = imread(fullfile(image_filenames(im_index).folder, image_filenames(im_index).name));
        im2 = imread(fullfile(image_filenames(im_index).folder, image_filenames(im_index + 1).name));
        
        % convert images to double arrays
        im1 = double(im1);
        im2 = double(im2);
        
        % mask images
        if io.image_masking
            im1 = im1 .* image_mask;
            im2 = im2 .* image_mask;
        end

        % flip images upside down
        im1 = flipud(im1);
        im2 = flipud(im2);
        image_mask = flipud(image_mask);
        
        if strcmp(io.image_type, 'experimental')
            % flip numbering of images if it is experimental. this is because in experimental data,
            % the order is gradient image followed by reference. the
            % displacement field will be flipped at the end to ensure
            % consistency.
            im_temp = im1;
            im1 = im2;
            im2 = im_temp;
        end
        %% load results from correlation processing for the current image pair

        if strcmp(tracking.initialization_method, 'correlation') || strcmp(id.diameter_estimation_method, 'correlation')
            % load vectors from cross-correlation
            vector_index = (image_pair_index-1)/(io.correlation_frame_step/2) + 1;
            if contains(vector_filenames(vector_index).name, 'corrplane')
                break;
            end
            correlation_results_current = load(fullfile(vector_filenames(vector_index).folder, vector_filenames(vector_index).name));
        end
        
        %% Identification and Sizing

        if strcmp(id.identification_method, 'standard')
            % ---------------------------------------------
            % standard dot identification and sizing
            % ---------------------------------------------
            
            %% Dot identification
            fprintf('running dot identification\n');

            % copy settings to data structure
%             particleIDprops = Data.ID;
            particleIDprops.v = id.intensity_threshold_current;
            particleIDprops.method = id.segmentation_method; %'dynamic';
            particleIDprops.contrast_ratio = 0;
            
            % perform dot identification for frame 1 (only for the first
            % image pair as all reference images are identical)
            if image_pair_index == 1
                [ID1_all{image_pair_index}.p_matrix,ID1_all{image_pair_index}.peaks,ID1_all{image_pair_index}.num_p]=particle_ID(im1,particleIDprops);
            else
                ID1_all{image_pair_index} = ID1_all{1};
            end
            
            % perform dot identification for frame 2
            [ID2_all{image_pair_index}.p_matrix,ID2_all{image_pair_index}.peaks,ID2_all{image_pair_index}.num_p]=particle_ID(im2,particleIDprops);
                
            %% Dot sizing
            fprintf('running dot sizing\n');

            % copy settings to data structure                        
            sizeprops.thresh = id.intensity_threshold_current;
            sizeprops.p_area = id.min_area; %d_p^2;
            sizeprops.sigma = 4;
            sizeprops.errors = double(sizing.default_iwc); % retain IWC estimate if gaussian fit fails
            sizeprops.method = sizing.centroid_subpixel_fit;
            
            % perform dot sizing for frame 1 (only for the first
            % image pair as all reference images are identical)
            if image_pair_index == 1
                [SIZE1_all{image_pair_index}.XYDiameter,SIZE1_all{image_pair_index}.mapsizeinfo,SIZE1_all{image_pair_index}.locxy]=particle_sizing(im1,ID1_all{image_pair_index}.p_matrix,...
                                    ID1_all{image_pair_index}.num_p,sizeprops);
                                
                % remove nan estimates from sizing results for frame 1
                [nan_r, ~] = find(isnan(SIZE1_all{image_pair_index}.XYDiameter));
                SIZE1_all{image_pair_index}.XYDiameter(nan_r, :) = [];
                SIZE1_all{image_pair_index}.mapsizeinfo(nan_r, :) = [];
                SIZE1_all{image_pair_index}.locxy(nan_r, :) = [];
                
                % extract x,y co-ordinates of identified dots in frame 2
                X1 = SIZE1_all{image_pair_index}.XYDiameter(:,1);
                Y1 = SIZE1_all{image_pair_index}.XYDiameter(:,2);
                Z1 = zeros(size(X1));
                
            else
                SIZE1_all{image_pair_index} = SIZE1_all{1};
            end
            
            % perform dot sizing for frame 2
            [SIZE2_all{image_pair_index}.XYDiameter,SIZE2_all{image_pair_index}.mapsizeinfo,SIZE2_all{image_pair_index}.locxy]=particle_sizing(im2,ID2_all{image_pair_index}.p_matrix,...
                                ID2_all{image_pair_index}.num_p,sizeprops);


            % remove nan estimates from sizing results for frame 2
            [nan_r, ~] = find(isnan(SIZE2_all{image_pair_index}.XYDiameter));
            SIZE2_all{image_pair_index}.XYDiameter(nan_r, :) = [];
            SIZE2_all{image_pair_index}.mapsizeinfo(nan_r, :) = [];
            SIZE2_all{image_pair_index}.locxy(nan_r, :) = [];

            
            % extract x,y co-ordinates of identified dots in frame 2
            X2 = SIZE2_all{image_pair_index}.XYDiameter(:,1);
            Y2 = SIZE2_all{image_pair_index}.XYDiameter(:,2);
            Z2 = zeros(size(X2));
        
        elseif contains(id.identification_method, 'apriori')
            % ------------------------------------------------------------
            % dot identification and sizing using prior information about
            % the dot locations
            % ------------------------------------------------------------

            %% calculate reference dot location from dot positions
            if image_pair_index == 1                
                
%                 [pos_ref_dots.x, pos_ref_dots.y] = calculate_reference_dot_locations_new(id.positions, experimental_parameters, id.camera_model, id.mapping_coefficients, id.order_z);
                pos_ref_dots.x = id.x_ref;
                pos_ref_dots.y = id.y_ref;
                num_dots_ref = numel(pos_ref_dots.x);
%                 if cropped_image
%                     pos_ref_dots.x = pos_ref_dots.x - crop_x;
%                     pos_ref_dots.x = NC_im - pos_ref_dots.x;
%                     pos_ref_dots.y = pos_ref_dots.y - crop_y;
%                 end
                
                [NR_im, NC_im] = size(im1);
                
                pos_ref_dots.y = NR_im - pos_ref_dots.y;
                
%                 % remove points outside FOV
%                 indices = pos_ref_dots.x < 1 | pos_ref_dots.x > NC_im-1 | pos_ref_dots.y < 1 | pos_ref_dots.y > NR_im-1;
%                 pos_ref_dots.x(indices) = [];
%                 pos_ref_dots.y(indices) = [];
%                 num_dots_ref = numel(pos_ref_dots.x); % - sum(indices);
                
                %if predicted locations lie in the masked region then ignore
                if io.image_masking
                    for dot_index = 1:num_dots_ref
                        %                 if image_mask(NR_im - round(pos_ref_dots.y(dot_index)), round(pos_ref_dots.x(dot_index))) == 0
                        if image_mask(round(pos_ref_dots.y(dot_index)), round(pos_ref_dots.x(dot_index))) == 0
                            pos_ref_dots.x(dot_index) = NaN;
                            pos_ref_dots.y(dot_index) = NaN;
                        end
                    end
                    nan_indices = isnan(pos_ref_dots.x) | isnan(pos_ref_dots.y);
                    pos_ref_dots.x(nan_indices) = [];
                    pos_ref_dots.y(nan_indices) = [];
                end
            end
            %% frame 1
            
            fprintf('im1\n');
            
            % estimate effective diameter from correlation
            if strcmp(id.diameter_estimation_method, 'correlation')
                d_p = estimate_effective_dot_diameter(correlation_results_current.X, correlation_results_current.Y, correlation_results_current.Di(:,:,end), pos_ref_dots.x, pos_ref_dots.y, id.dot_diameter);
            else
                d_p = id.dot_diameter*ones(size(pos_ref_dots.x));
            end
            
            % set minimum area for a set of pixels to be considered a dot
            min_area = 0.5 * median(d_p)^2;
            
            % identify and size dots using their known locations on the
            % target (only for first image pair, as all reference images
            % are identical)
            if image_pair_index == 1
                [SIZE1_all{image_pair_index}.XYDiameter,SIZE1_all{image_pair_index}.peaks, SIZE1_all{image_pair_index}.mapsizeinfo,SIZE1_all{image_pair_index}.locxy, SIZE1_all{image_pair_index}.mapint]=combined_ID_size_apriori_07(im1,pos_ref_dots.x, pos_ref_dots.y, d_p+2, sizing.centroid_subpixel_fit, sizing.default_iwc, id.min_area, id.W_area, id.W_intensity, id.W_distance);                

                % extract co-ordinates
                X1 = SIZE1_all{image_pair_index}.XYDiameter(:,1);
                Y1 = SIZE1_all{image_pair_index}.XYDiameter(:,2);
                Z1 = zeros(size(X1));
            else
                SIZE1_all{image_pair_index} = SIZE1_all{1};
            end
            
            %% frame 2
            
            fprintf('im2\n');
            
            % estimate probable locations of dots on frame 2 based on their
            % position in frame 1 and the displacement field from
            % correlation or light ray displacements
            if strcmp(id.initialization_method, 'rays')
                [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, id.initialization_method, pos_ref_rays_1.x, pos_ref_rays_1.y, pos_ref_rays_2.x - pos_ref_rays_1.x, pos_ref_rays_2.y - pos_ref_rays_1.y);
            elseif strcmp(id.initialization_method, 'correlation')
                [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, id.initialization_method, correlation_results_current.X, correlation_results_current.Y, correlation_results_current.U, correlation_results_current.V);
            else
                X2_est = X1;
                Y2_est = Y1;
            end
                       
            % identify and size dots using their known locations on the target
            [SIZE2_all{image_pair_index}.XYDiameter,SIZE2_all{image_pair_index}.peaks,SIZE2_all{image_pair_index}.mapsizeinfo,SIZE2_all{image_pair_index}.locxy, SIZE2_all{image_pair_index}.mapint]=combined_ID_size_apriori_07(im2, X2_est, Y2_est, d_p+2, sizing.centroid_subpixel_fit, sizing.default_iwc, id.min_area, id.W_area, id.W_intensity, id.W_distance);
            
            % extract co-ordinates
            X2 = SIZE2_all{image_pair_index}.XYDiameter(:,1);
            Y2 = SIZE2_all{image_pair_index}.XYDiameter(:,2);
            Z2 = zeros(size(X2));
        else
            fprintf('unknown input\n');
            keyboard;
        end
        
        %% Tracking
        
        fprintf('running dot tracking\n');
        % if a guess for the dot positions in the second frame is required
        % and does not already exist, then calculae it
        if ~strcmp(tracking.initialization_method, 'none')
            if ~exist('X2_est', 'var') && ~exist('Y2_est', 'var')            
                if strcmp(image_type, 'original')
                    [X2_est, Y2_est] = calculate_predicted_dot_positions(X1, Y1, tracking.initialization_method, correlation_results_current.X, correlation_results_current.Y, correlation_results_current.U, correlation_results_current.V);
                elseif strcmp(image_type, 'deform')
                    X2_est = X1;
                    Y2_est = Y1;
                end
            end
        % if the guess is not required then set it to be the positions of 
        % the dots in the first frame
        else
            X2_est = X1;
            Y2_est = Y1;
        end
        
        % expected z co-ordinates of dots in the second frame (zero)
        Z2_est = zeros(size(X1));

        % extract dot diameters from frame 1
        d1 = SIZE1_all{image_pair_index}.XYDiameter(:,3);
        % extract dot diameters from frame 2
        d2 = SIZE2_all{image_pair_index}.XYDiameter(:,3);
        
        % extract dot intensities from frame 1
        I1 = SIZE1_all{image_pair_index}.XYDiameter(:,4);
        % extract dot intensities from frame 2
        I2 = SIZE2_all{image_pair_index}.XYDiameter(:,4);

        % track the particles in the image pair using the 3D weighted
        % nearest neighbor tracking method
        [tracks_all{image_pair_index}]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,...
            Z1,Z2,Z2_est,d1,d2,I1,I2, [tracking.distance_weight, tracking.size_weight, tracking.intensity_weight], tracking.search_radius);

        %% correct sub-pixel estimate of displacement using a correlation based estimate
        if tracking.perform_correlation_correction
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
                % extract current track
                track_current = tracks_all{image_pair_index}(track_index, :);
                [U(track_index), V(track_index)] = cross_correlate_dots_06(im1, im2, SIZE1_all{image_pair_index}, SIZE2_all{image_pair_index}, track_current, tracking.correlation_correction.algorithm, tracking.correlation_correction.subpixel_fit, tracking.correlation_correction.zero_mean, tracking.correlation_correction.min_sub);            
            end
            % append results to track
            tracks_all{image_pair_index} = [tracks_all{image_pair_index}, U, V];
        end
        
        %% flip tracking results to account for reordering reference and gradient images
        
        if strcmp(io.image_type, 'experimental')
            tracks_temp = tracks_all{image_pair_index};
            % swap x positions
            tracks_temp(:, [1, 2]) = tracks_all{image_pair_index}(:, [2, 1]);
            % swap y positions
            tracks_temp(:, [3, 4]) = tracks_all{image_pair_index}(:, [4, 3]);
            % swap z positions
            tracks_temp(:, [5, 6]) = tracks_all{image_pair_index}(:, [6, 5]);
            % swap diameters
            tracks_temp(:, [7, 8]) = tracks_all{image_pair_index}(:, [8, 7]);
            % swap intensities
            tracks_temp(:, [9, 10]) = tracks_all{image_pair_index}(:, [10, 9]);
            % swap particle ids
            tracks_temp(:, [11, 12]) = tracks_all{image_pair_index}(:, [12, 11]);
            
            if tracking.perform_correlation_correction
                % change sign on displacements from correlation correction
                tracks_temp(:, 14) = -tracks_all{image_pair_index}(:, 14);
                tracks_temp(:, 15) = -tracks_all{image_pair_index}(:, 15);
            end
            
            % save to original data structure
            tracks_all{image_pair_index} = tracks_temp;
            
            % flip id results
            ID1_temp = ID1_all{image_pair_index};
            ID1_all{image_pair_index} = ID2_all{image_pair_index};
            ID2_all{image_pair_index} = ID1_temp;
            
            % flip sizing results
            SIZE1_temp = SIZE1_all{image_pair_index};
            SIZE1_all{image_pair_index} = SIZE2_all{image_pair_index};
            SIZE2_all{image_pair_index} = SIZE1_temp;            
            
        end
        
        %% save results to file
        
        save_filename = ['file_' num2str(im_index, '%04d') '.mat'];
        % save identification results
        id1 = ID1_all{image_pair_index};
        id2 = ID2_all{image_pair_index};
        save(fullfile(id_save_directory, save_filename), 'id1', 'id2');
        
        % save sizing results
        size1 = SIZE1_all{image_pair_index};
        size2 = SIZE2_all{image_pair_index};
        save(fullfile(size_save_directory, save_filename), 'size1', 'size2');
        
        % save tracking
        tracks = tracks_all{image_pair_index};
        save(fullfile(track_save_directory, save_filename), 'tracks');
        
    end
    
end
