clear
close all
clc

restoredefaultpath;

addpath ../dot-tracking-package/
dbstop if error

logical_string = {'False', 'True'};

%% experiment settings
% date of experiment
test_date = '2019-03-22';
% electrode gap (mm)
gap = 2;
% view (front, side, top)
camera_view = 'front-view';
% top level project directory
project_directory = '/scratch/shannon/c/aether/Projects/plasma-induced-flow/';% body placed in the tunnel
% top directory containing images for this case
top_image_directory = fullfile(project_directory, 'analysis', 'data', 'dbd', test_date, [num2str(gap) 'mm'], camera_view, 'bos'); 
% top directory to contain results for this case
top_results_directory = fullfile(project_directory, 'analysis', 'results', 'dbd', test_date, [num2str(gap) 'mm'], camera_view, 'bos');
% pulse voltage (V)
pulse_voltage = 400;
% pulse width (ns)
pulse_width = 40;

% dot diameter in the image plane (pix.)
dot_diameter_image = 4;
% dot diameter in the object plane (mm)
dot_diameter_object = 0.5;

% dot spacing in the image plane (pix.)
dot_spacing_image = 1;
% magnification (um/pix.)
magnification = dot_diameter_object * 1e3/dot_diameter_image; %33;

% size of the camera sensor (um)
experimental_parameters.camera_design.pixel_pitch = 20.8; %20;
% number of pixels on the camera
experimental_parameters.camera_design.x_pixel_number = 512; %1000;
experimental_parameters.camera_design.y_pixel_number = 352; %309;
% camera angles (deg.)
experimental_parameters.camera_design.x_camera_angle = 0;
experimental_parameters.camera_design.y_camera_angle = 0;
% distance between camera and the dot target (um)
% (not required if camera model is soloff)
experimental_parameters.lens_design.object_distance = NaN; 
% non-dimensional magnification
experimental_parameters.lens_design.magnification = experimental_parameters.camera_design.pixel_pitch/magnification; %6.45/28.63; %1.3;
% spacing between dots (um)
experimental_parameters.bos_pattern.dot_spacing = dot_spacing_image * magnification; %dot_spacing * magnification; % 2*150; %2*250; %180;

%% Processing settings

% ------------------------------------------
% I/O settings
% ------------------------------------------
io = struct;

% image index starting
io.imfstart = 1;
% time interval between correlating image pairs
io.imcstep = 1;
% number of image pairs to skip during processing
io.imfstep = 2;
% last image to be read
io.imfend = 1;

% perform image masking (true/false)
io.image_masking = true;
% image_mask_filepath = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-05/0_25mm/allpsi/test1/processing/single-ref/staticmask.tif';
io.image_mask_filepath = fullfile(top_image_directory, 'mask.tif');

% were the images cropped prior to processing? (true/false)
io.cropped_image = false;
% image type (synthetic or experimental)
io.image_type = 'experimental';

% directory containing results from correlation analysis (if hybrid
% tracking or diameter estimation from peak sizing is to be done)
io.correlation_results.read_directory = '';
% type of multi pass scheme used in the correlation
io.correlation_results.multi_pass_scheme = 'deform';
% pass from which results are to be used for initializing the tracking
io.correlation_results.pass_number = 2;
% number of images (not pairs) to skip during correlation
io.correlation_results.frame_step = 2;
% this is the base name of the files that contain the results to be analyzed
io.correlation_results.basename =  ['BOS*pass' num2str(io.correlation_results.pass_number) '_'];
% ------------------------------------------
% Identification and Sizing settings
% ------------------------------------------
id = struct;
% image segmentation method (for standard id)
id.segmentation_method = 'blob';
% intensity threshold for standard identification
id.intensity_threshold_current = 500; 

% position estimation method for the reference image ('apriori' or
% 'standard')
id.identification_method = 'apriori';

% subpixel fit ('tpg', 'lsg')
sizing.centroid_subpixel_fit = 'iwc'; %'clsg';
% default to iwc if the gaussian fit fails? (true/false)
sizing.default_iwc = false;

% ------------------------------------------
% Apriori position estimation settings
% ------------------------------------------

% approximate dot diameter [pix.]
id.dot_diameter = dot_diameter_image;
% minimum area for a set of pixels to be considered a dot [pix.^2]
id.min_area = 0.5 * id.dot_diameter^2;
% camera model to use for projecting known dot positions from object space
% to image space if apriori identification is to be performed
% ('thin-lens', 'soloff')
id.camera_model = 'soloff';
% order of the z polynomial mapping function for soloff (if application)
id.order_z = 1;
% directory containing the camera model
id.camera_model_directory = fullfile(top_image_directory, 'calibration');
% method to estimate approximate effective diameter of a dot in the image
% ('correlation', 'none');
id.diameter_estimation_method = 'none';
% weights for separating true peak from false peaks in dynamic segmentation
id.W_area = 1;
id.W_intensity = 1;
id.W_distance = 1;
% offset on the camera sensor (pix.)
id.starting_index_x = 0;
id.starting_index_y = 1;
% initialization method for hybrid tracking ('none', 'correlation')
id.initialization_method = 'none';

% ------------------------------------------
% Tracking settings
% ------------------------------------------
tracking = struct;
% search radius for nearest neighbor search [pix.]
tracking.search_radius = 5;
% weights for the nearest neighbor algorithm
tracking.distance_weight = 1; %1/3;
tracking.size_weight = 0; %1/3;
tracking.intensity_weight = 0; %1/3;
% type of image on which tracking is to be performed (original or deformed)
tracking.image_type = 'original';
% initialization method for hybrid tracking ('none', 'correlation')
tracking.initialization_method = 'none';

% ------------------------------------------
% Correlation Correction settings
% ------------------------------------------

% perform correlation_correction? (true/false)
tracking.perform_correlation_correction = true;
% correlation correction algorithm ('dcc', 'scc', 'rpc')
tracking.correlation_correction.algorithm = 'dcc';
% zero mean correlation windows? (true/false)
tracking.correlation_correction.zero_mean = 0;
% perform minimim subtraction? (true/false)
tracking.correlation_correction.min_sub = 1;
% subpixel fit for the correlation plane ('tpg', 'lsg')
tracking.correlation_correction.subpixel_fit = 'lsg';

%% load sample job file

% this is the filepath containing the sample parameter file
sample_job_filepath = '/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/sample-job-files/';

% this is the sample parameter filename 
sample_job_filename = ['sample_job_ptv_' sizing.centroid_subpixel_fit '.mat'];

% load a sample parameter file
sample_bos_params = load([sample_job_filepath sample_job_filename]);
Data = sample_bos_params.Data;
Data.par = '0';

%% process images

% obtain list of folders in this directory
[shot_folders, num_shots] = get_directory_listing(fullfile(top_image_directory, [num2str(pulse_voltage) 'V_' num2str(pulse_width) 'ns']), 'shot*');

for shot_index = 1 %:num_shots
    %% create directory to store results
    % display progress to user
    fprintf('run number: %d\n', shot_index);

    % top read directory for this case
    current_image_directory = fullfile(shot_folders(shot_index).folder, shot_folders(shot_index).name, 'processing', 'average-ref');
    current_results_directory = fullfile(top_results_directory, [num2str(pulse_voltage) 'V_' num2str(pulse_width) 'ns'], shot_folders(shot_index).name, 'ptv');

    % directory containing images
    io.image_directory = current_image_directory;

    % name for this processing case
    case_name = ['id=' id.identification_method '_size=' sizing.centroid_subpixel_fit '_corr=' tracking.correlation_correction.subpixel_fit];

    % create directory to save results
    io.results_save_directory = fullfile(current_results_directory, case_name);
    if ~exist(io.results_save_directory, 'dir')
        mkdir(io.results_save_directory);
    end

    % run tracking
    run_dot_tracking_v2(io, id, sizing, tracking, experimental_parameters);
end
