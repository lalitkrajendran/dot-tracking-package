clear
close all
clc

dbstop if error

%% set image generation parameters

% geometric dotsize (um)
d_geo = 0;
% diffraction diameter (pix.)
d_diff = 3.00;
% pixel pitch (um)
pixel_pitch = 10;
% focal length (mm)
focal_length = 105;
% f-number
f_number = 11;
% presence of overlapping dots
dot_overlap = false;
% cropping
crop_x = 256;
crop_y = 256;

% seeding density
seeding_density_current = 20;
% noise levels
noise_std_array = [0];

imzeros = '5';

% pixel locations of the origin marker
x0 = 536;
y0 = 485;

% update pixel location to account for the first pixel being 1 and not 0
starting_index_x = -0.5; %-0.5;
starting_index_y = -0.5; %-0.5;

% no overlap case
top_image_directory = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/jhu/overlap=False/general/512x512/rk4/zobj1000-zmin245-zmax255-nx256-ny256-nz64/';
num_image_pairs = 10;
cropped_image = false;

%% set processing parameters
% --------------------------------
% TRACKING
% --------------------------------

% fraction by which the intensity threshold is to be reduced each iteration
threshold_relaxation_parameter = 0.1;

% flag to turn on image smoothing
image_smoothing = false;

% position estimation method for the reference image ('apriori' or
% 'standard')
% method = 'standard';
% method = 'iterative';
method = 'apriori-dots-size-07'; %'apriori-dots-size-04'
camera_model = 'thin-lens';
order_z = 2;
% method used to create residual image in iterative ID
dot_removal_method = 2;

% method to estimate dot diameter on the image for use in identification
% can be 1) true (synthetic images), 2) correlation (from correlation peak 
% size) or 3) none
diameter_estimation_method = 'true';

% tolerance to remove ghost points based on known dot positions
tolerance = 1;

% subpixel fit
subpixel_fit = 'lsg';
% flag to select if the iwc estimate is to be taken as the default if the
% gaussian fit fails
default_iwc = false;

% correlation correction for final displacement
correlation_correction = true;
logical_string = {'False'; 'True'};
correlation_correction_algorithm = 'dcc-v5-min-sub';
zero_mean = 0;
min_sub = 1;
correlation_correction_subpixel_fit = 'tpg';

% search radius for nearest neighbor search (pix.)
search_radius = 5;
% weights for the nearest neighbor algorithm
distance_weight = 1; %1/3;
size_weight = 0; %1/3;
intensity_weight = 0; %1/3;

% weights for separating true peak from false peaks in dynamic segmentation
W_area = 1;
W_intensity = 1;
W_distance = 1;

% --------------------------------
% CORRELATION 
% --------------------------------
window = false;
% initialization using correlation results
correlation_initialization = false; %true;
% intialization method for tracking ('rays' or 'correlation' or 'none')
initialization_method = 'none'; %'rays';
multi_pass_scheme = 'deform';
window_overlap_percentage = '50';

correlation_results_filename = 'BOS_pass2_*.mat';
% case name for reading correlation results
case_name = 'rk4-scc-2pass-deform-32-16-overlap50-tpg-validated';
% pass number for correlation results
pass_number = 2;

%% load sample job file

% this is the filepath containing the sample parameter file
sample_job_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/sample-job-files/';

% this is the sample parameter filename 
sample_job_filename = ['sample_job_ptv_' subpixel_fit '.mat'];

% load a sample parameter file
sample_bos_params = load([sample_job_filepath sample_job_filename]);
Data = sample_bos_params.Data;

%% adjust settings for job file
Data.par = '0';
Data.imzeros = '5';

% --------------------------
% ID settings
% --------------------------
Data.ID.method = 'dynamic';
% Data.ID.v = intensity_threshold;
Data.ID.run = 1;
Data.ID.contrast_ratio = 0;
particleIDprops = Data.ID;

% --------------------------
% Size settings
% --------------------------
Data.Size.run      = 1;
% Data.Size.thresh   = intensity_threshold;
Data.Size.method   = subpixel_fit; %sizing_method_current;
Data.Size.size_threshold_fraction = 0.5;
% Data.Size.p_area   = d_diff^2 * Data.Size.size_threshold_fraction; %3;
Data.Size.p_area   = 3;
Data.Size.sigma    = 4;
Data.Size.errors   = 1;
sizeprops = Data.Size;
min_area = 8;

%% process images

% time_array = [0,8]; %, 8];
time_array = [8, 12, 12, 14, 14]; %, 8];
zoff_array = [0.75, 0.5, 0.75, 0.25, 0.75];

num_snapshots = length(time_array);

for snapshot_index = 1:num_snapshots
time = time_array(snapshot_index);
zoff = zoff_array(snapshot_index);
fprintf('time: %.2f, zoff: %.2f pi \n', time, zoff);    %% modify read/write settings for job file

for noise_std_index = 1:length(noise_std_array)
    noise_std_current = noise_std_array(noise_std_index);
    fprintf('noise: %d\n', noise_std_current);
    
    % directory containing images for the current case
    current_image_directory = fullfile(top_image_directory, ['f=' num2str(focal_length) 'mm_f' num2str(f_number, '%02d') '_l_p=' num2str(pixel_pitch) 'um'], ['dg=' num2str(d_geo) 'um_dp=' num2str(d_diff, '%.2f') 'pix'], ['seeding=' num2str(seeding_density_current, '%.2f')], ['t=' num2str(time, '%.2f') '_zoff=' num2str(zoff, '%.2f') 'pi'], ['noise' num2str(noise_std_current, '%02d')]);
    
    if contains(method, 'apriori')        
        results_folder_name_1 = fullfile(current_image_directory, 'processing', 'results', 'ptv', method, camera_model, ['diameter_estimation=' diameter_estimation_method]);
        str_split = strsplit(method, '-');
        method_version = str2double(str_split{end});        
        if method_version >= 4
            results_folder_name_1 = fullfile(results_folder_name_1, ['id-weights=' num2str(W_area, '%.2f'), '_' num2str(W_intensity, '%.2f'), '_' num2str(W_distance, '%.2f')]);
        end
    else
        results_folder_name_1 = fullfile(current_image_directory, 'processing', 'results', 'ptv', method);
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

    %% load data for initializing tracking
    filenames = [];
    if strcmp(initialization_method, 'correlation')
        % ------------------------
        % load correlation results
        % ------------------------
        correlation_results_filepath = fullfile(current_image_directory, 'processing', case_name, 'results', 'vectors');
        % this populates all the files in the directory
        filenames = dir(fullfile(correlation_results_filepath, correlation_results_filename));
%         data = load(fullfile(correlation_results_filepath, correlation_results_filename));
    elseif strcmp(initialization_method, 'rays')
        % ------------------------
        % calculate reference displacements from light ray positions
        % ------------------------
        ray_displacements_all.x = [];
        ray_displacements_all.y = [];
        for image_pair_index = 1:num_image_pairs
            light_ray_data_filepath = fullfile(current_image_directory, '1');
            [pos_1, pos_2, dir_1, dir_2] = load_lightray_data_02(light_ray_data_filepath);
            [ray_displacements, ray_deflections] = calculate_lightray_deflections(pos_1, pos_2, dir_1, dir_2);
            ray_displacements_all.x = [ray_displacements_all.x; ray_displacements.x];
            ray_displacements_all.y = [ray_displacements_all.y; ray_displacements.y];
        end
        % find average ray deflection along x and y (microns)
        del_x = nanmean(ray_displacements_all.x);
        del_y = nanmean(ray_displacements_all.y);

        % convert to pixel units
        del_x_pixels = del_x/pixel_pitch;
        del_y_pixels = del_y/pixel_pitch;

        del_x_std_pixels = nanstd(ray_displacements_all.x)/pixel_pitch;

%         % copy the values to the variables already in the code
%         U_ref(delta_x_index, noise_std_index) = del_x_pixels;
    end
    
    %% process image pairs

    ID1_temp = cell(1,num_image_pairs);
    ID2_temp = cell(1,num_image_pairs);
    SIZE1_temp = cell(1,num_image_pairs);
    SIZE2_temp = cell(1,num_image_pairs);
    
    % calcualte intensity threshold for identification
    intensity_threshold_current = 3000; %noise_std_current * 0.001 * 65535;
    Data.ID.v = intensity_threshold_current;
    Data.Size.thresh   = intensity_threshold_current;
    Data.Size.p_area = d_diff^2;
    
    sizeprops = Data.Size;
    sizeprops.p_area = d_diff^2;
    
    particleIDprops = Data.ID;
    
    pos_ref = cell(1,4);

    for image_pair_index = 1:num_image_pairs
        fprintf('image pair index: %d\n', image_pair_index);
        % load parameters used to generate image
        parameters = load(fullfile(current_image_directory,  num2str(image_pair_index), 'parameters.mat'));
        % load positions of dots on the dot pattern in the world co-ordinates
        positions = load(fullfile(current_image_directory,  num2str(image_pair_index), 'positions.mat'));        
        % load images
        imdirec = fullfile(current_image_directory, num2str(image_pair_index));
        im1 = imread(fullfile(imdirec, 'bos_pattern_image_1.tif'));
        im2 = imread(fullfile(imdirec, 'bos_pattern_image_2.tif'));

        %% pre-process image
%         im1 = pre_process_image(im1);
%         im2 = pre_process_image(im2);

        %% smooth image
        if image_smoothing
%             im1 = imfilter(im1, 1/9*ones(3));
%             im2 = imfilter(im2, 1/9*ones(3));
            im1 = medfilt2(im1, [3, 3]);
            im2 = medfilt2(im2, [3, 3]);
        end
        
        %% calculate dot positions from light rays

        % load light ray positions
        [pos_1, pos_2, dir_1, dir_2] = load_lightray_data_02(fullfile(current_image_directory, num2str(image_pair_index)));            
        [x_pos_1, y_pos_1] = pos_to_xy_pix_03(pos_1, parameters, starting_index_x, starting_index_y);
        [x_pos_2, y_pos_2] = pos_to_xy_pix_03(pos_2, parameters, starting_index_x, starting_index_y);

        % calculate centroid from light rays
        num_lightrays_per_dot = parameters.bos_pattern.particle_number_per_grid_point*parameters.bos_pattern.lightray_number_per_particle ;
%         [pos_ref_rays_1.x, pos_ref_rays_1.y] = calculate_centroids_from_lightrays_04(x_pos_1, y_pos_1);
%         [pos_ref_rays_2.x, pos_ref_rays_2.y] = calculate_centroids_from_lightrays_04(x_pos_2, y_pos_2);

%         pos_ref_rays_1 = sort_positions(pos_ref_rays_1);
%         pos_ref_rays_2 = sort_positions(pos_ref_rays_2);

%         [pos_ref_rays_1.x, pos_ref_rays_1.y, indices_1] = calculate_centroids_from_lightrays_05(x_pos_1, y_pos_1, d_diff);
%         [pos_ref_rays_2.x, pos_ref_rays_2.y, indices_2] = calculate_centroids_from_lightrays_05(x_pos_2, y_pos_2, d_diff);

        [pos_ref_rays_1.x, pos_ref_rays_1.y] = calculate_centroids_from_lightrays_06(x_pos_1, y_pos_1, num_lightrays_per_dot);
        [pos_ref_rays_2.x, pos_ref_rays_2.y] = calculate_centroids_from_lightrays_06(x_pos_2, y_pos_2, num_lightrays_per_dot);

        if cropped_image
            pos_ref_rays_1.x = pos_ref_rays_1.x - crop_x;
            pos_ref_rays_1.y = pos_ref_rays_1.y - crop_y;

            pos_ref_rays_2.x = pos_ref_rays_2.x - crop_x;
            pos_ref_rays_2.y = pos_ref_rays_2.y - crop_y;
        end
        
        pos_ref_rays_1.x = parameters.camera_design.x_pixel_number - pos_ref_rays_1.x;
        pos_ref_rays_2.x = parameters.camera_design.y_pixel_number - pos_ref_rays_2.x;

        % remove nan values
        nan_indices = isnan(pos_ref_rays_1.x) | isnan(pos_ref_rays_1.y) | isnan(pos_ref_rays_2.x) | isnan(pos_ref_rays_2.y);
        pos_ref_rays_1.x(nan_indices) = [];
        pos_ref_rays_1.y(nan_indices) = [];
        pos_ref_rays_2.x(nan_indices) = [];
        pos_ref_rays_2.y(nan_indices) = [];
        
        %% load correlation results for initializing tracking
        
        if strcmp(initialization_method, 'correlation')
            % load vectors from cross-correlation
            data = load(fullfile(correlation_results_filepath, filenames(image_pair_index).name));
            [NR, NC] = size(data.U);

            data.Y = parameters.camera_design.y_pixel_number - data.Y;
            data.V = -data.V;
            
%             % account for sign convention used in the code
%             data.V(:) = -data.V(:);
%             data.U(:) = -data.U;
        end
        
        %% calculate reference dot location from dot positions
        if strcmp(camera_model, 'soloff')
            % load camera mapping coefficients
            mapping_coefficients = load(fullfile(calibration_directory, ['camera_model_type=' num2str(order_z) '.mat']));
        else
            mapping_coefficients = [];
        end
        [pos_ref_dots.x, pos_ref_dots.y] = calculate_reference_dot_locations_new(positions, parameters, camera_model, mapping_coefficients, order_z, starting_index_x, starting_index_y);
        pos_ref_dots = sort_positions(pos_ref_dots);

        if cropped_image
            pos_ref_dots.x = pos_ref_dots.x - crop_x;
            pos_ref_dots.y = pos_ref_dots.y - crop_y;
        end
        
        pos_ref_dots.x = parameters.camera_design.x_pixel_number - pos_ref_dots.x;

        % align the first co-ordinates
        pos_ref_dots.x = pos_ref_dots.x + min(pos_ref_rays_1.x) - min(pos_ref_dots.x); %(pos_ref_rays.x(1) - pos_ref_dots.x(1));
        pos_ref_dots.y = pos_ref_dots.y + min(pos_ref_rays_1.y) - min(pos_ref_dots.y); %(pos_ref_rays.y(1) - pos_ref_dots.y(1));

        indices = pos_ref_dots.x < 0 | pos_ref_dots.x > parameters.camera_design.x_pixel_number | pos_ref_dots.y < 0 | pos_ref_dots.y > parameters.camera_design.y_pixel_number;
%             pos_ref_dots.x(indices) = [];
%             pos_ref_dots.y(indices) = [];
        num_dots_ref = numel(pos_ref_dots.x) - sum(indices);

        %% Identification and Sizing        
        
        if strcmp(method, 'iterative')
            [X1, Y1] = identify_dots_iterative(im1, Data, parameters, positions, 'ref', dot_removal_method, threshold_relaxation_parameter, num_dots_ref, d_diff/2);
            [X2, Y2] = identify_dots_iterative(im2, Data, parameters, positions, 'grad', dot_removal_method, threshold_relaxation_parameter, num_dots_ref, d_diff/2);

            %% remove ghost points
    %             [X1, Y1] = remove_ghost_points(X1, Y1, pos_ref_dots.x, pos_ref_dots.y, tolerance);
            [X1, Y1] = remove_ghost_points(X1, Y1, pos_ref_rays_1.x, pos_ref_rays_1.y, d_diff/2);
    %             X1 = X1_iter;
    %             Y1 = Y1_iter;
        elseif strcmp(method, 'standard')
            %% Dot identification

            fprintf('running dot identification\n');

            [ID1_temp{image_pair_index}.p_matrix,ID1_temp{image_pair_index}.peaks,ID1_temp{image_pair_index}.num_p]=particle_ID(im1,particleIDprops);
            [ID2_temp{image_pair_index}.p_matrix,ID2_temp{image_pair_index}.peaks,ID2_temp{image_pair_index}.num_p]=particle_ID(im2,particleIDprops);
            
            %% Dot sizing
            fprintf('running dot sizing\n');

            % pos_prana = cell(1,num_sizing_methods);
            pos_prana = [];

            [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy]=particle_sizing(im1,ID1_temp{image_pair_index}.p_matrix,...
                                ID1_temp{image_pair_index}.num_p,sizeprops);
            [SIZE2_temp{image_pair_index}.XYDiameter,SIZE2_temp{image_pair_index}.mapsizeinfo,SIZE2_temp{image_pair_index}.locxy]=particle_sizing(im2,ID2_temp{image_pair_index}.p_matrix,...
                                ID2_temp{image_pair_index}.num_p,sizeprops);

            X1 = SIZE1_temp{image_pair_index}.XYDiameter(:,1);
            Y1 = SIZE1_temp{image_pair_index}.XYDiameter(:,2);

            X2 = SIZE2_temp{image_pair_index}.XYDiameter(:,1);
            Y2 = SIZE2_temp{image_pair_index}.XYDiameter(:,2);

        else
            if contains(method, 'apriori')
                %% image 1
                
                fprintf('im1\n');
%                 % identify dots using their known locations on the target
%                 if contains(method, 'rays')
%                     [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy]=combined_ID_size_apriori_02(im1,pos_ref_rays_1.x, pos_ref_rays_1.y, d_diff, subpixel_fit);
%                 elseif contains(method, 'dots')
%                     [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy, SIZE1_temp{image_pair_index}.mapint]=combined_ID_size_apriori_03(im1,pos_ref_dots.x, pos_ref_dots.y, d_diff+2, subpixel_fit, min_area);
%                 end

                % identify dots using their known locations on the target
                if contains(method, 'rays')
                    [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy]=combined_ID_size_apriori_02(im1,pos_ref_rays_1.x, pos_ref_rays_1.y, d_diff, subpixel_fit);
                elseif contains(method, 'dots')
%                     [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy, SIZE1_temp{image_pair_index}.mapint]=combined_ID_size_apriori_03(im1,pos_ref_dots.x, pos_ref_dots.y, d_diff+2, subpixel_fit, min_area);
%                     [SIZE1_temp{image_pair_index}.XYDiameter, peaks, SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy, SIZE1_temp{image_pair_index}.mapint]=combined_ID_size_apriori_04(im1,pos_ref_dots.x, pos_ref_dots.y, d_diff+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                    [SIZE1_temp{image_pair_index}.XYDiameter,SIZE1_temp{image_pair_index}.peaks, SIZE1_temp{image_pair_index}.mapsizeinfo,SIZE1_temp{image_pair_index}.locxy, SIZE1_temp{image_pair_index}.mapint]=combined_ID_size_apriori_07(im1,pos_ref_dots.x, pos_ref_dots.y, d_diff+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                end

                X1 = SIZE1_temp{image_pair_index}.XYDiameter(:,1);
                Y1 = SIZE1_temp{image_pair_index}.XYDiameter(:,2);
                %% image 2

                if correlation_initialization
                    if strcmp(initialization_method, 'rays')                    
                        [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, pos_ref_rays_1.x, pos_ref_rays_1.y, pos_ref_rays_2.x - pos_ref_rays_1.x, pos_ref_rays_2.y - pos_ref_rays_1.y);                
                    elseif strcmp(initialization_method, 'correlation')
                        [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, data.X, data.Y, data.U, data.V);                
                    end
                else
                    X2_est = X1;
                    Y2_est = Y1;
                end
                
                fprintf('im2\n');
                
%                 [SIZE2_temp{image_pair_index}.XYDiameter,SIZE2_temp{image_pair_index}.mapsizeinfo,SIZE2_temp{image_pair_index}.locxy, SIZE2_temp{image_pair_index}.mapint]=combined_ID_size_apriori_03(im2, X2_est, Y2_est, d_diff+2, subpixel_fit, min_area);
%                 [SIZE2_temp{image_pair_index}.XYDiameter, peaks, SIZE2_temp{image_pair_index}.mapsizeinfo,SIZE2_temp{image_pair_index}.locxy, SIZE2_temp{image_pair_index}.mapint]=combined_ID_size_apriori_04(im2, X2_est, Y2_est, d_diff+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                [SIZE2_temp{image_pair_index}.XYDiameter,SIZE2_temp{image_pair_index}.peaks,SIZE2_temp{image_pair_index}.mapsizeinfo,SIZE2_temp{image_pair_index}.locxy, SIZE2_temp{image_pair_index}.mapint]=combined_ID_size_apriori_07(im2, X2_est, Y2_est, d_diff+2, subpixel_fit, default_iwc, min_area, W_area, W_intensity, W_distance);
                
                X2 = SIZE2_temp{image_pair_index}.XYDiameter(:,1);
                Y2 = SIZE2_temp{image_pair_index}.XYDiameter(:,2);

            else
                fprintf('unknown input\n');
                keyboard;
            end
        end
        ID1_temp{image_pair_index}.num_dots_ref = num_dots_ref;

        % set z co-ordinates to zero
        Z1 = zeros(size(X1));
        Z2 = zeros(size(X2));
        
        % calculate position error
%         [pos_err_x, pos_err_y] = calculate_position_error(X1, Y1, pos_ref_rays_1.x, pos_ref_rays_1.y, tolerance, true, d_diff);
        %% Tracking
        fprintf('running dot tracking\n');

        if correlation_initialization && ~exist('X2_est') && ~exist('Y2_est')            
%             [X2_est, Y2_est] = calculate_predicted_dot_positions(X1, Y1, initialization_method, data.X, data.Y, data.U, data.V, []);                
            if strcmp(initialization_method, 'rays')                    
                [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, pos_ref_rays_1.x, pos_ref_rays_1.y, pos_ref_rays_2.x - pos_ref_rays_1.x, pos_ref_rays_2.y - pos_ref_rays_1.y);                
            elseif strcmp(initialization_method, 'correlation')
                [X2_est, Y2_est] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, data.X, data.Y, data.U, data.V);                
            end       
        end

        Z2_est = Z1;
        
        d1 = SIZE1_temp{image_pair_index}.XYDiameter(:,3);
        d2 = SIZE2_temp{image_pair_index}.XYDiameter(:,3);
        I1 = SIZE1_temp{image_pair_index}.XYDiameter(:,4);
        I2 = SIZE2_temp{image_pair_index}.XYDiameter(:,4);
        
        weights = [distance_weight, size_weight, intensity_weight]; %[1, 0, 0];

        s_radius = search_radius;
        
%         plot_predicted_dot_positions(X1, Y1, X2, Y2, X2_est, Y2_est);
        % run weighted nearest neighbors
        %track the particles in the image pair using the 3D weighted
        %nearest neighbor tracking method
        [tracks_temp{image_pair_index}]=weighted_nearest_neighbor3D(X1,X2,X2_est,Y1,Y2,Y2_est,...
            Z1,Z2,Z2_est,d1,d2,I1,I2,weights,s_radius);
        
%         % calculate position error
%         [pos_err_x, pos_err_y] = calculate_position_error(tracks_temp{image_pair_index}(:,1), tracks_temp{image_pair_index}(:,3), pos_ref_rays_1.x, pos_ref_rays_1.y, tolerance, false, d_diff);
%         [pos_err_x, pos_err_y] = calculate_position_error(tracks_temp{image_pair_index}(:,2), tracks_temp{image_pair_index}(:,4), pos_ref_rays_2.x, pos_ref_rays_2.y, tolerance, false, d_diff);
        
        % plot position error histogram
        
        %% correct sub-pixel estimate of displacement using a correlation based estimate
        if correlation_correction
            % get final sub-pixel displacement estimate by cross-correlation intensity
            % maps of the dots
            num_tracks = size(tracks_temp{image_pair_index}, 1);
            U = zeros(num_tracks,1);
            V = zeros(num_tracks,1);
            for track_index = 1:num_tracks
    %             fprintf('track number: %d\n', track_index);
                track_current = tracks_temp{image_pair_index}(track_index, :);
%                 [U(track_index), V(track_index)] = cross_correlate_dots(im1, im2, track_current(1), track_current(2), track_current(3), track_current(4), track_current(7), track_current(8), correlation_correction_algorithm, subpixel_fit);            
%                 [U(track_index), V(track_index)] = cross_correlate_dots_02(im1, im2, SIZE1_temp{image_pair_index}, SIZE2_temp{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit);            
%                 [U(track_index), V(track_index)] = cross_correlate_dots_03(SIZE1_temp{image_pair_index}, SIZE2_temp{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit);            
                [U(track_index), V(track_index)] = cross_correlate_dots_05(SIZE1_temp{image_pair_index}, SIZE2_temp{image_pair_index}, track_current, correlation_correction_algorithm, correlation_correction_subpixel_fit, zero_mean, min_sub);            
            end
            % append results to track
            tracks_temp{image_pair_index} = [tracks_temp{image_pair_index}, U, V];
        end
        
        fprintf('saving results\n');
        % ----------------------
        % save SIZE results
        % ----------------------
        
        save_filepath = fullfile(current_results_directory, 'size');
        if ~exist(save_filepath, 'dir')
            mkdir(save_filepath);
        end
        save_filename = fullfile(save_filepath, ['ptv_bos_Size_' num2str(image_pair_index*2 - 1, ['%0' imzeros 'd']), '.mat']);
        Size = SIZE1_temp{image_pair_index};
        save(save_filename, 'Size');

        save_filename = fullfile(save_filepath, ['ptv_bos_Size_' num2str(image_pair_index*2, ['%0' imzeros 'd']), '.mat']);
        Size = SIZE2_temp{image_pair_index};
        save(save_filename, 'Size');

        % ----------------------
        % save tracking results
        % ----------------------
        save_filepath = fullfile(current_results_directory, 'track');
        if ~exist(save_filepath, 'dir')
            mkdir(save_filepath);
        end
        save_filename = fullfile(save_filepath, ['ptv_bos_Track_' num2str(image_pair_index*2 - 1, ['%0' imzeros 'd']), '.mat']);
        tracks = tracks_temp{image_pair_index};
        save(save_filename, 'tracks');

    end    
end
end


    