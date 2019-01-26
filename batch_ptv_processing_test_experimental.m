clear
close all
clc

% % add prana to path
% addpath('/home/barracuda/a/lrajendr/Software/prana/');
% addpath /home/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/crlb/
% addpath /home/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/error-analysis-codes/ptv-codes/

dbstop if error

%% case settings

% top directory containing images
% top_image_directory = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-05/0_25mm/allpsi/';
top_image_directory = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-11/0_15mm/allpsi/';

skip_image_pairs = 1;
correlation_frame_step = 2;

% image masking
image_masking = true;
% image_mask_filepath = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-05/0_25mm/allpsi/test1/processing/single-ref/staticmask.tif';
image_mask_filepath = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-11/0_15mm/allpsi/test1/processing/single-ref/staticmask.tif';
% position estimation method for the reference image ('apriori' or
% 'standard')
method = 'apriori-dots-size-07';
% method = 'standard';
% camera model
% camera_model = 'thin-lens';
camera_model = 'soloff';
order_z = 1;
% directory containing the camera model
camera_model_directory = top_image_directory;

% method to estimate approximate effective diameter of a dot in the image
% diameter_estimation_method = 'none';
diameter_estimation_method = 'correlation';

% type of image on which tracking is to be performed (original or deformed)
image_type = 'original';

% subpixel fit
subpixel_fit = 'tpg';
% flag to select if the iwc estimate is to be taken as the default if the
% gaussian fit fails
default_iwc = false;

% weights for separating true peak from false peaks in dynamic segmentation
W_area = 1;
W_intensity = 1;
W_distance = 1;

% list of particle identification methods
IDmethods = {'blob','dynamic','combined'}; 

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

% initialization_method = 'correlation';
initialization_method = 'none';

multi_pass_scheme = 'deform';

correlation_correction = true;
correlation_correction_algorithm = 'dcc-v6-min-sub';
zero_mean = 1;
min_sub = 0;
correlation_correction_subpixel_fit = 'lsg';
logical_string = {'False', 'True'};

% pass from which results are to be used for initializing the tracking
pass_number = 2;

% this is the base name of the files that contain the results to be analyzed
results_basename =  ['BOS*pass' num2str(pass_number) '_'];

cropped_image = false;

x_align = 6;
y_align = 4;

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

    % display progress to user
    fprintf('test number: %d\n', test_index);

    % top read directory for this case
    current_image_directory = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'single-ref');

    % directory containing processing results for ptv images    
    if contains(method, 'apriori')        
        results_folder_name_1 = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'results', 'ptv', method, camera_model, ['diameter_estimation=' diameter_estimation_method]);
        if contains(method, '04')
            results_folder_name_1 = fullfile(results_folder_name_1, ['id-weights=' num2str(W_area, '%.2f'), '_' num2str(W_intensity, '%.2f'), '_' num2str(W_distance, '%.2f')]);
        end
    else
        results_folder_name_1 = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'results', 'ptv', method);
    end

    results_folder_name_2 = fullfile(subpixel_fit, ['default_iwc=' logical_string{default_iwc + 1}], ['s_radius=' num2str(search_radius)], ...
            ['weights=' num2str(distance_weight, '%.2f'), '_' num2str(size_weight, '%.2f'), '_' num2str(intensity_weight, '%.2f')]);
        
    if strcmp(initialization_method, 'correlation')
        results_folder_name_2 = fullfile(results_folder_name_2, ['initialization-' initialization_method '-' multi_pass_scheme]);
    else
        results_folder_name_2 = fullfile(results_folder_name_2, ['initialization-' initialization_method]);
    end
    
    if correlation_correction
        results_folder_name_3 = fullfile(['correlation_correction=' logical_string{correlation_correction+1}], ['correction_algorithm=' correlation_correction_algorithm '-' correlation_correction_subpixel_fit]);
    else
        results_folder_name_3 = fullfile(['correlation_correction=' logical_string{correlation_correction+1}]);
    end

    current_results_directory = fullfile(results_folder_name_1, results_folder_name_2, results_folder_name_3);
    
    if ~exist(current_results_directory, 'dir')
        mkdir(current_results_directory);
    end

    % copy the parameters
    Data = sample_bos_params.Data;

    % get directory listing for images
    [image_filenames,num_image_files] = get_directory_listing(current_image_directory, 'im*.tif' );
    
    % modify number of images to read
    num_image_pairs = num_image_files/2; 
    if strcmp(initialization_method, 'correlation')
        vectors_directory = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'vectors-single-reference-deform');
        % get directory listing for vectors
        [vector_filenames,~] = get_directory_listing(vectors_directory, [results_basename '*.mat']);
    end
    
    if strcmp(diameter_estimation_method, 'correlation')
        diameters_directory = fullfile(top_image_directory, test_folders(test_index).name, 'processing', 'vectors-single-reference-deform');
        % get directory listing for vectors
        [diameter_filenames,~] = get_directory_listing(diameters_directory, [results_basename '*.mat']);
        
    end
    
    if image_masking
%         mask_img = imread('/home/shannon/c/aether/Projects/BOS/spark-induced-flow/analysis/data/10-30-2017/staticmask.tif');
        image_mask = imread(image_mask_filepath);
        % convert mask to 0 and 1
        image_mask(image_mask > 0) = 1;
        image_mask = double(image_mask);
    end
    
    ID1_all = cell(1,num_image_pairs);
    ID2_all = cell(1,num_image_pairs);
    SIZE1_all = cell(1,num_image_pairs);
    SIZE2_all = cell(1,num_image_pairs);
    tracks_all = cell(1,num_image_pairs);
    Data_all = cell(1,num_image_pairs);
    
    %% calculate reference dot location from dot positions
    if contains(method, 'apriori')
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
            data = load(fullfile(camera_model_directory, 'calibration_data.mat'));
            parameters.bos_pattern.X_Min = min(data.calibration_data.x_world_full{1}) * 1e3;
            parameters.bos_pattern.Y_Min = min(data.calibration_data.y_world_full{1}) * 1e3;
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
%         pos_ref_dots = sort_positions(pos_ref_dots);
        if cropped_image
            pos_ref_dots.x = pos_ref_dots.x - crop_x;
            pos_ref_dots.x = NC_im - pos_ref_dots.x;
            pos_ref_dots.y = pos_ref_dots.y - crop_y;
        end
%         pos_ref_dots.x = pos_ref_dots.x - min(pos_ref_dots.x);
%         pos_ref_dots.y = pos_ref_dots.y - min(pos_ref_dots.y);
        
        % align the first co-ordinates
        % pos_ref_dots_all.x = pos_ref_dots_all.x + (pos_ref_rays_all.x(1) - pos_ref_dots_all.x(1));
        % pos_ref_dots_all.y = pos_ref_dots_all.y + (pos_ref_rays_all.y(1) - pos_ref_dots_all.y(1));
%         pos_ref_dots.x = pos_ref_dots.x + x_align - min(pos_ref_dots.x); %(pos_ref_rays.x(1) - pos_ref_dots.x(1));
%         pos_ref_dots.y = pos_ref_dots.y + y_align - min(pos_ref_dots.y); %(pos_ref_rays.y(1) - pos_ref_dots.y(1));

        % load a sample image to get its dimensions
        im = imread(fullfile(current_image_directory, image_filenames(1).name));
        [NR_im, NC_im] = size(im);
        
%         pos_ref_dots.y = NR_im - pos_ref_dots.y;
        
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
        a = 1;
    elseif strcmp(method, 'calibration')
        % load camera mapping coefficients
        mapping_coefficients = load(fullfile(camera_model_directory, ['camera_model_type=' num2str(order_z) '.mat']));
        % load data used for calibration
        calibration_data = load(fullfile(camera_model_directory, 'calibration_data.mat'));
    end
    

    % calculate intensity threshold for identification
    intensity_threshold_current = 5000; %100; %noise_std_current * 0.001 * 65535;
    Data.ID.v = intensity_threshold_current;
    Data.Size.thresh   = intensity_threshold_current;
    Data.Size.p_area = min_area; %d_p^2;
    Data.Size.sigma = 4;
    Data.Size.errors = 1; % retain IWC estimate if gaussian fit fails
    sizeprops = Data.Size;
    particleIDprops = Data.ID;
    particleIDprops.method = 'blob'; %'dynamic';
    particleIDprops.contrast_ratio = 0;
    for image_pair_index = 1:skip_image_pairs:num_image_pairs %2:num_image_pairs        
        Data_all{image_pair_index} = Data;
        fprintf('Image Pair Number: %d\n', image_pair_index);
        %% calculate dot positions from prana

        % load images
        im1 = imread(fullfile(current_image_directory, image_filenames((image_pair_index-1)*2 + 1).name));
        im2 = imread(fullfile(current_image_directory, image_filenames((image_pair_index-1)*2 + 2).name));
        
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
        % displacement field will be flipped at the end
        im_temp = im1;
        im1 = im2;
        im2 = im_temp;
        %% load correlation results
        
        if strcmp(initialization_method, 'correlation') 
            % load vectors from cross-correlation
            vector_index = (image_pair_index-1)/(correlation_frame_step/2) + 1;
            if contains(vector_filenames(vector_index).name, 'corrplane')
                break;
            end
            data = load(fullfile(vectors_directory, vector_filenames(vector_index).name));
%             data.Y = size(im1,1) - data.Y;
%             data.V = -data.V;
        end
        
        if strcmp(diameter_estimation_method, 'correlation')
            % load vectors from cross-correlation
            diameter_index = (image_pair_index-1)/(correlation_frame_step/2) + 1;
            if contains(diameter_filenames(diameter_index).name, 'corrplane')
                break;
            end
            diameter_data = load(fullfile(diameters_directory, diameter_filenames(diameter_index).name));                            
        end
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
        tracks_temp{image_pair_index} = tracks_temp;
        
        %% save tracks
        fprintf('Saving the tracks\n');
        save_filepath = fullfile(current_results_directory, 'track');
        if ~exist(save_filepath, 'dir')
            mkdir(save_filepath);
        end

        vector_index = (image_pair_index-1)/(correlation_frame_step/2) + 1;
        % save tracking results
%         save_filename = fullfile(save_filepath, ['ptv_bos_Track_' num2str(image_pair_index*2 - 1, ['%0' Data_temp{image_pair_index}.imzeros 'd']), '.mat']);
        save_filename = fullfile(save_filepath, ['ptv_bos_Track_' num2str(image_pair_index*2 - 1, '%05d'), '.mat']);
        tracks = tracks_all{image_pair_index};
        save(save_filename, 'tracks');
    end
end
    