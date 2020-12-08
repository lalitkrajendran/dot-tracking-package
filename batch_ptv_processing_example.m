clear
close all
clc

restoredefaultpath;

dbstop if error

logical_string = {'False', 'True'};

% ==========================
%% experiment settings
% ==========================
% dot diameter in the image plane (pix.)
dot_diameter = 7;
% dot spacing in the image plane (pix.)
dot_spacing = 8;
% magnification (um/pix.)
magnification = 10.5; %33;

% size of the camera sensor (um)
experimental_parameters.camera_design.pixel_pitch = 20.8; %20;
% number of pixels on the camera
experimental_parameters.camera_design.x_pixel_number = 1024; %1000;
experimental_parameters.camera_design.y_pixel_number = 1024; %309;
% camera angles (deg.)
experimental_parameters.camera_design.x_camera_angle = 0;
experimental_parameters.camera_design.y_camera_angle = 0;
% distance between camera and the dot target (um)
experimental_parameters.lens_design.object_distance = 10.5 * 2.54 * 10e3;
% non-dimensional magnification
experimental_parameters.lens_design.magnification = experimental_parameters.camera_design.pixel_pitch/magnification; %6.45/28.63; %1.3;
% spacing between dots (um)
experimental_parameters.bos_pattern.dot_spacing = 84; %dot_spacing * magnification; % 2*150; %2*250; %180;

% ==========================
%% Processing settings
% ==========================
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
io.image_mask_filepath = './sample-data/mask.tif';

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
% display intermediate progress to user? (true/false)
io.display_intermediate_progress = false;

% ------------------------------------------
% Standard Identification settings
% ------------------------------------------
id = struct;
% image segmentation method (for standard id)
id.segmentation_method = 'blob';
% intensity threshold for standard identification
id.intensity_threshold_current = 500; 

% position estimation method for the reference image ('apriori' or
% 'standard')
id.identification_method = 'apriori';

% ------------------------------------------
% Apriori position estimation settings
% ------------------------------------------
% approximate dot diameter [pix.]
id.dot_diameter = dot_diameter;
% minimum area for a set of pixels to be considered a dot [pix.^2]
id.min_area = 0.5 * id.dot_diameter^2;
% camera model to use for projecting known dot positions from object space
% to image space if apriori identification is to be performed
% ('thin-lens', 'soloff')
id.camera_model = 'soloff';
% order of the z polynomial mapping function for soloff (if application)
id.order_z = 1;
% % directory containing the camera model
% id.camera_model_directory = fullfile(top_image_directory, 'calibration');
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
% Sizing settings for both standard/apriori
% ------------------------------------------
% subpixel fit ('tpg', 'lsg')
sizing.centroid_subpixel_fit = 'lsg'; 
% default to iwc if the gaussian fit fails? (true/false)
sizing.default_iwc = false;

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

% ------------------------------------------
% Validation settings
% ------------------------------------------
% perform validation? (True/False)
tracking.perform_validation = true;
% perform displacement thresholding? (True/False)
tracking.validation.perform_displacement_thresholding = true;
% maximum displacement threshold (pix.)
tracking.validation.max_displacement_threshold = 3;
% perform UOD? (True/False)
tracking.validation.perform_uod = true;
% replace vectors in UOD? (True/False)
tracking.validation.replace_vectors = true;
% UOD residual threshold (2-4, lower is stricter)
tracking.validation.uod_residual_threshold = 2;
% minimum allowable threshold [pix.]
tracking.validation.uod_epsilon = 0.1;

% ==========================
%% process images
% ==========================
% top read directory for this case
current_image_directory = './sample-data/images/'; 
current_results_directory = './sample-results/';

% directory containing images
io.image_directory = current_image_directory;

% directory containing calibration results
id.camera_model_directory = fullfile(current_results_directory, 'calibration');

% name for this processing case
case_name = ['id=' id.identification_method '_size=' sizing.centroid_subpixel_fit '_corr=' tracking.correlation_correction.subpixel_fit];
% create directory to save results
io.results_save_directory = fullfile(current_results_directory, case_name);
if ~exist(io.results_save_directory, 'dir')
    mkdir(io.results_save_directory);
end

% run tracking
run_dot_tracking_v5(io, id, sizing, tracking, experimental_parameters);

