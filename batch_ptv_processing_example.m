clear
close all
clc

dbstop if error

logical_string = {'False', 'True'};
%% Case settings

% ------------------------------------------
% I/O settings
% ------------------------------------------

% top directory containing images
top_image_directory = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-11/0_15mm/allpsi/';
% number of image pairs to skip during processing
image_pair_increment = 1;
% number of images (not pairs) to skip during correlation
correlation_frame_step = 2;

% perform image masking (true/false)
image_masking = true;
% image_mask_filepath = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-05/0_25mm/allpsi/test1/processing/single-ref/staticmask.tif';
image_mask_filepath = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-11/0_15mm/allpsi/test1/processing/single-ref/staticmask.tif';

% were the images cropped prior to processing? (true/false)
cropped_image = false;

% ------------------------------------------
% Identification and Sizing settings
% ------------------------------------------

% dot identification method for standard identification
% ('blob','dynamic','combined')
ID_method = 'blob'; 
% intensity threshold for standard identification
intensity_threshold_current = 5000; 

% position estimation method for the reference image ('apriori' or
% 'standard')
method = 'apriori';

% subpixel fit ('tpg', 'lsg')
subpixel_fit = 'tpg';
% default to iwc if the gaussian fit fails? (true/false)
default_iwc = false;

% ------------------------------------------
% Apriori position estimation settings
% ------------------------------------------

% camera model to use for projecting known dot positions from object space
% to image space if apriori identification is to be performed
% ('thin-lens', 'soloff')
camera_model = 'soloff';
% order of the z polynomial mapping function for soloff (if application)
order_z = 1;
% directory containing the camera model
camera_model_directory = top_image_directory;
% method to estimate approximate effective diameter of a dot in the image
% ('correlation', 'none');
diameter_estimation_method = 'correlation';
% weights for separating true peak from false peaks in dynamic segmentation
W_area = 1;
W_intensity = 1;
W_distance = 1;

% ------------------------------------------
% Tracking settings
% ------------------------------------------

% search radius for nearest neighbor search [pix.]
search_radius = 5;
% approximate dot diameter [pix.]
d_p_approx = 10;
% minimum area for a set of pixels to be considered a dot [pix.^2]
min_area = 0.5 * d_p_approx^2;
% weights for the nearest neighbor algorithm
distance_weight = 1; %1/3;
size_weight = 0; %1/3;
intensity_weight = 0; %1/3;
% type of image on which tracking is to be performed (original or deformed)
image_type = 'original';
% initialization method for hybrid tracking ('none', 'correlation')
initialization_method = 'none';

% ------------------------------------------
% Hybrid Tracking settings
% ------------------------------------------

% type of multi pass scheme used in the correlation
multi_pass_scheme = 'deform';
% pass from which results are to be used for initializing the tracking
pass_number = 2;
% this is the base name of the files that contain the results to be analyzed
results_basename =  ['BOS*pass' num2str(pass_number) '_'];

% ------------------------------------------
% Correlation Correction settings
% ------------------------------------------

% perform correlation_correction? (true/false)
correlation_correction = true;
% correlation correction algorithm ('dcc', 'scc', 'rpc')
correlation_correction_algorithm = 'dcc';
% zero mean correlation windows? (true/false)
zero_mean = 0;
% perform minimim subtraction? (true/false)
min_sub = 1;
% subpixel fit for the correlation plane ('tpg', 'lsg')
correlation_correction_subpixel_fit = 'lsg';

%% load sample job file

% this is the filepath containing the sample parameter file
sample_job_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/sample-job-files/';

% this is the sample parameter filename 
sample_job_filename = ['sample_job_ptv_' subpixel_fit '.mat'];

% load a sample parameter file
sample_bos_params = load([sample_job_filepath sample_job_filename]);
Data = sample_bos_params.Data;
Data.par = '0';

%% process images

% get number of tests available
[test_folders, num_test_folders] = get_directory_listing(top_image_directory, 'test*');

for test_index = 1:num_test_folders
    %% create directory to store results
    % display progress to user
    fprintf('test number: %d\n', test_index);

    % top read directory for this case
    current_image_directory = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'single-ref');
    
    % first part of the results directory name 
    if contains(method, 'apriori')        
        results_folder_name_1 = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'results', 'ptv', method, camera_model, ['diameter_estimation=' diameter_estimation_method]);
        if contains(method, '04')
            results_folder_name_1 = fullfile(results_folder_name_1, ['id-weights=' num2str(W_area, '%.2f'), '_' num2str(W_intensity, '%.2f'), '_' num2str(W_distance, '%.2f')]);
        end
    else
        results_folder_name_1 = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'results', 'ptv', method);
    end

    % second part of the results directory name    
    results_folder_name_2 = fullfile(subpixel_fit, ['default_iwc=' logical_string{default_iwc + 1}], ['s_radius=' num2str(search_radius)], ...
            ['weights=' num2str(distance_weight, '%.2f'), '_' num2str(size_weight, '%.2f'), '_' num2str(intensity_weight, '%.2f')]);
        
    if strcmp(initialization_method, 'correlation')
        results_folder_name_2 = fullfile(results_folder_name_2, ['initialization-' initialization_method '-' multi_pass_scheme]);
    else
        results_folder_name_2 = fullfile(results_folder_name_2, ['initialization-' initialization_method]);
    end

    % third part of the results directory name    
    if correlation_correction
        results_folder_name_3 = fullfile(['correlation_correction=' logical_string{correlation_correction+1}], ['correction_algorithm=' correlation_correction_algorithm '-' correlation_correction_subpixel_fit]);
    else
        results_folder_name_3 = fullfile(['correlation_correction=' logical_string{correlation_correction+1}]);
    end

    % full name of the results directory
    current_results_directory = fullfile(results_folder_name_1, results_folder_name_2, results_folder_name_3);
    % create directory if it does not exist
    if ~exist(current_results_directory, 'dir')
        mkdir(current_results_directory);
    end

    %% get directory listings for images, vectors and diameters

    % get directory listing for images
    [image_filenames,num_image_files] = get_directory_listing(current_image_directory, 'im*.tif' );    
    % calculate number of image pairs to read
    num_image_pairs = num_image_files/2; 
    
    % get directory listing for correlation vectors results
    if strcmp(initialization_method, 'correlation') || strcmp(diameter_estimation_method, 'correlation')
        vectors_directory = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'vectors-single-reference-deform');
        % get directory listing for vectors
        [vector_filenames,~] = get_directory_listing(vectors_directory, [results_basename '*.mat']);
    end
        
    if image_masking
        % load image mask
        image_mask = imread(image_mask_filepath);
        % convert mask to 0 and 1
        image_mask(image_mask > 0) = 1;
        % convert mask to double
        image_mask = double(image_mask);
    end
    
    %% load and process image pairs

    % initialize data structures to hold results
    ID1_all = cell(1,num_image_pairs);
    ID2_all = cell(1,num_image_pairs);
    SIZE1_all = cell(1,num_image_pairs);
    SIZE2_all = cell(1,num_image_pairs);
    tracks_all = cell(1,num_image_pairs);

    for image_pair_index = 1:image_pair_increment:num_image_pairs %2:num_image_pairs        
        % display progress to user
        fprintf('Image Pair Number: %d\n', image_pair_index);
        
        %% load image pairs and perform masking

        % load images
        im1 = imread(fullfile(current_image_directory, image_filenames((image_pair_index-1)*2 + 1).name));
        im2 = imread(fullfile(current_image_directory, image_filenames((image_pair_index-1)*2 + 2).name));
        
        % convert images to double arrays
        im1 = double(im1);
        im2 = double(im2);
        
        % mask images
        if image_masking
            im1 = im1 .* image_mask;
            im2 = im2 .* image_mask;
        end

        % flip images upside down
        im1 = flipud(im1);
        im2 = flipud(im2);
        
        % flip numbering of images. this is because in experimental data,
        % the order is gradient image followed by reference. the
        % displacement field will be flipped at the end to ensure
        % consistency.
        im_temp = im1;
        im1 = im2;
        im2 = im_temp;
        
        %% load results from correlation processing for the current image pair

        if strcmp(initialization_method, 'correlation') || strcmp(diameter_estimation_method, 'correlation')
            % load vectors from cross-correlation
            vector_index = (image_pair_index-1)/(correlation_frame_step/2) + 1;
            if contains(vector_filenames(vector_index).name, 'corrplane')
                break;
            end
            correlation_results_current = load(fullfile(vectors_directory, vector_filenames(vector_index).name));
        end
        
        %% Identification and Sizing

        if strcmp(method, 'standard')
            % ---------------------------------------------
            % standard dot identification and sizing
            % ---------------------------------------------
            
            %% Dot identification
            fprintf('running dot identification\n');

            % copy settings to data structure
            particleIDprops = Data.ID;
            particleIDprops.v = intensity_threshold_current;
            particleIDprops.method = ID_method; %'dynamic';
            particleIDprops.contrast_ratio = 0;
            
            % perform dot identification for frame 1
            [ID1_all{image_pair_index}.p_matrix,ID1_all{image_pair_index}.peaks,ID1_all{image_pair_index}.num_p]=particle_ID(im1,particleIDprops);
            % perform dot identification for frame 1
            [ID2_all{image_pair_index}.p_matrix,ID2_all{image_pair_index}.peaks,ID2_all{image_pair_index}.num_p]=particle_ID(im2,particleIDprops);
                
            %% Dot sizing
            fprintf('running dot sizing\n');

            % copy settings to data structure            
            sizeprops = Data.Size;
            sizeprops.thresh = intensity_threshold_current;
            sizeprops.p_area = min_area; %d_p^2;
            sizeprops.sigma = 4;
            sizeprops.errors = double(default_iwc); % retain IWC estimate if gaussian fit fails

            % perform dot sizing for frame 1
            [SIZE1_all{image_pair_index}.XYDiameter,SIZE1_all{image_pair_index}.mapsizeinfo,SIZE1_all{image_pair_index}.locxy]=particle_sizing(im1,ID1_all{image_pair_index}.p_matrix,...
                                ID1_all{image_pair_index}.num_p,sizeprops);
            % perform dot sizing for frame 2
            [SIZE2_all{image_pair_index}.XYDiameter,SIZE2_all{image_pair_index}.mapsizeinfo,SIZE2_all{image_pair_index}.locxy]=particle_sizing(im2,ID2_all{image_pair_index}.p_matrix,...
                                ID2_all{image_pair_index}.num_p,sizeprops);

            % remove nan estimates from sizing resluts for frame 1
            [nan_r, ~] = find(isnan(SIZE1_all{image_pair_index}.XYDiameter));
            SIZE1_all{image_pair_index}.XYDiameter(nan_r, :) = [];
            SIZE1_all{image_pair_index}.mapsizeinfo(nan_r, :) = [];
            SIZE1_all{image_pair_index}.locxy(nan_r, :) = [];

            % remove nan estimates from sizing resluts for frame 1
            [nan_r, ~] = find(isnan(SIZE2_all{image_pair_index}.XYDiameter));
            SIZE2_all{image_pair_index}.XYDiameter(nan_r, :) = [];
            SIZE2_all{image_pair_index}.mapsizeinfo(nan_r, :) = [];
            SIZE2_all{image_pair_index}.locxy(nan_r, :) = [];

            % extract x,y co-ordinates of identified dots in frame 1
            X1 = SIZE1_all{image_pair_index}.XYDiameter(:,1);
            Y1 = SIZE1_all{image_pair_index}.XYDiameter(:,2);
            Z1 = zeros(size(X1));
            
            % extract x,y co-ordinates of identified dots in frame 2
            X2 = SIZE2_all{image_pair_index}.XYDiameter(:,1);
            Y2 = SIZE2_all{image_pair_index}.XYDiameter(:,2);
            Z2 = zeros(size(X2));

        elseif contains(method, 'apriori')
            % ------------------------------------------------------------
            % dot identification and sizing using prior information about
            % the dot locations
            % ------------------------------------------------------------

            %% calculate reference dot location from dot positions
            if image_pair_index == 1
                % camera parameters (all units in um)
                parameters.camera_design.pixel_pitch = 6.45; %20;
                parameters.camera_design.x_pixel_number = 1392; %1000;
                parameters.camera_design.y_pixel_number = 1040; %309;
                parameters.camera_design.x_camera_angle = 0;
                parameters.camera_design.y_camera_angle = 0;
                % lens parameters
                parameters.lens_design.object_distance = 22e4; %6 * 2.54 * 10e3;
                parameters.lens_design.magnification = 6.45/28.63; %6.45/28.63; %1.3;
                
                % bos pattern parameters
                if strcmp(camera_model, 'soloff')
                    calibration_results = load(fullfile(camera_model_directory, 'calibration_data.mat'));
                    parameters.bos_pattern.X_Min = min(calibration_results.calibration_data.x_world_full{1}) * 1e3;
                    parameters.bos_pattern.Y_Min = min(calibration_results.calibration_data.y_world_full{1}) * 1e3;
                else
                    parameters.bos_pattern.X_Min = 0; %-1.5e3; %0;
                    parameters.bos_pattern.Y_Min = 0; %-3e3; %0;
                end
                % spacing between dots (um)
                parameters.bos_pattern.dot_spacing = 2*150; %2*250; %180;
                
                field_of_view = parameters.camera_design.x_pixel_number * parameters.camera_design.pixel_pitch / parameters.lens_design.magnification;
                parameters.bos_pattern.dot_number = round(field_of_view/parameters.bos_pattern.dot_spacing) * 2;
                
                x_array = parameters.bos_pattern.X_Min + parameters.bos_pattern.dot_spacing * (0:parameters.bos_pattern.dot_number);
                y_array = parameters.bos_pattern.Y_Min + parameters.bos_pattern.dot_spacing * (0:parameters.bos_pattern.dot_number);
                
                [positions.x, positions.y] = meshgrid(x_array, y_array);
                
                if strcmp(camera_model, 'soloff')
                    % load camera mapping coefficients
                    mapping_coefficients = load(fullfile(camera_model_directory, ['camera_model_type=' num2str(order_z) '.mat']));
                else
                    mapping_coefficients = [];
                end
                
                [pos_ref_dots.x, pos_ref_dots.y] = calculate_reference_dot_locations_new(positions, parameters, camera_model, mapping_coefficients, order_z);
                
                if cropped_image
                    pos_ref_dots.x = pos_ref_dots.x - crop_x;
                    pos_ref_dots.x = NC_im - pos_ref_dots.x;
                    pos_ref_dots.y = pos_ref_dots.y - crop_y;
                end
                
                % load a sample image to get its dimensions
                im = imread(fullfile(current_image_directory, image_filenames(1).name));
                [NR_im, NC_im] = size(im);
                
                pos_ref_dots.y = NR_im - pos_ref_dots.y;
                
                % remove points outside FOV
                indices = pos_ref_dots.x < 1 | pos_ref_dots.x > NC_im-1 | pos_ref_dots.y < 1 | pos_ref_dots.y > NR_im-1;
                pos_ref_dots.x(indices) = [];
                pos_ref_dots.y(indices) = [];
                num_dots_ref = numel(pos_ref_dots.x); % - sum(indices);
                
                % if predicted locations lie in the masked region then ignore
                if image_masking
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
            if strcmp(diameter_estimation_method, 'correlation')
                d_p = estimate_effective_dot_diameter(correlation_results_current.X, correlation_results_current.Y, correlation_results_current.Di(:,:,end), pos_ref_dots.x, pos_ref_dots.y, d_p_approx);
            else
                d_p = d_p_approx*ones(size(pos_ref_dots.x));
            end
            
            % set minimum area for a set of pixels to be considered a dot
            min_area = 0.5 * median(d_p)^2;
            
            % identify and size dots using their known locations on the target
            [SIZE1_all{image_pair_index}.XYDiameter,SIZE1_all{image_pair_index}.peaks, SIZE1_all{image_pair_index}.mapsizeinfo,SIZE1_all{image_pair_index}.locxy, SIZE1_all{image_pair_index}.mapint]=combined_ID_size_apriori_07(im1,pos_ref_dots.x, pos_ref_dots.y, d_p+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);

            % extract co-ordinates
            X1 = SIZE1_all{image_pair_index}.XYDiameter(:,1);
            Y1 = SIZE1_all{image_pair_index}.XYDiameter(:,2);
            Z1 = zeros(size(X1));
            %% frame 2
            
            fprintf('im2\n');
            
            % estimate probable locations of dots on frame 2 based on their
            % position in frame 1 and the displacement field from
            % correlation or light ray displacements
            if strcmp(initialization_method, 'rays')
                [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, pos_ref_rays_1.x, pos_ref_rays_1.y, pos_ref_rays_2.x - pos_ref_rays_1.x, pos_ref_rays_2.y - pos_ref_rays_1.y);
            elseif strcmp(initialization_method, 'correlation')
                [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, correlation_results_current.X, correlation_results_current.Y, correlation_results_current.U, correlation_results_current.V);
            else
                X2_est = X1;
                Y2_est = Y1;
            end
            
            % estimate effective diameter from correlation at the guessed
            % dot locations in frame 2
            if strcmp(diameter_estimation_method, 'correlation')
                d_p = estimate_effective_dot_diameter(correlation_results_current.X, correlation_results_current.Y, correlation_results_current.Di(:,:,end), X2_est, Y2_est, d_p_approx);
            else
                d_p = d_p_approx*ones(size(pos_ref_dots.x));
            end
            
            % identify and size dots using their known locations on the target
            [SIZE2_all{image_pair_index}.XYDiameter,SIZE2_all{image_pair_index}.peaks,SIZE2_all{image_pair_index}.mapsizeinfo,SIZE2_all{image_pair_index}.locxy, SIZE2_all{image_pair_index}.mapint]=combined_ID_size_apriori_07(im2, X2_est, Y2_est, d_p+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
            
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
        if ~strcmp(initialization_method, 'none')
            if ~exist('X2_est', 'var') && ~exist('Y2_est', 'var')            
                if strcmp(image_type, 'original')
                    [X2_est, Y2_est] = calculate_predicted_dot_positions(X1, Y1, initialization_method, correlation_results_current.X, correlation_results_current.Y, correlation_results_current.U, correlation_results_current.V);
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
            Z1,Z2,Z2_est,d1,d2,I1,I2, [distance_weight, size_weight, intensity_weight], search_radius);

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
                % extract current track
                track_current = tracks_all{image_pair_index}(track_index, :);
                [U(track_index), V(track_index)] = cross_correlate_dots_06(im1, im2, SIZE1_all{image_pair_index}, SIZE2_all{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit, zero_mean, min_sub);            
            end
            % append results to track
            tracks_all{image_pair_index} = [tracks_all{image_pair_index}, U, V];
        end
        
        %% flip tracking results to account for reordering reference and gradient images
        tracks_temp = tracks_all{image_pair_index};
        % swap x positions
        tracks_temp(:, [1, 2]) = tracks_temp(:, [2, 1]);
        % swap y positions
        tracks_temp(:, [3, 4]) = tracks_temp(:, [4, 3]);
        % swap z positions
        tracks_temp(:, [5, 6]) = tracks_temp(:, [6, 5]);
        % swap diameters
        tracks_temp(:, [7, 8]) = tracks_temp(:, [8, 7]);
        % swap intensities
        tracks_temp(:, [9, 10]) = tracks_temp(:, [10, 9]);
        % swap particle ids
        tracks_temp(:, [11, 12]) = tracks_temp(:, [12, 11]);
        
        if correlation_correction
            % change sign on displacements from correlation correction
            tracks_temp(:, 14) = -tracks_temp(:, 14);
            tracks_temp(:, 15) = -tracks_temp(:, 15);
        end
        
        % save to original data structure
        tracks_all{image_pair_index} = tracks_temp;
        
        %% save tracks
        fprintf('Saving the tracks\n');
        save_filepath = fullfile(current_results_directory, 'track');
        if ~exist(save_filepath, 'dir')
            mkdir(save_filepath);
        end

        vector_index = (image_pair_index-1)/(correlation_frame_step/2) + 1;
        % save tracking results
        save_filename = fullfile(save_filepath, ['ptv_bos_Track_' num2str(image_pair_index*2 - 1, '%05d'), '.mat']);
        tracks = tracks_all{image_pair_index};
        save(save_filename, 'tracks');
    end
end
    