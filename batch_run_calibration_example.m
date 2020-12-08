% example script to calculate mapping function coefficients from reference image

clear
close all
clc

restoredefaultpath;

% ==========================
%% calibration settings
% ==========================
% dot size (mm)
dot_size = 0.042;
% diameter of a dot
dot_diameter = dot_size * 1e-3;
% spacing between dots
dot_spacing = 2 * dot_size * 1e-3;
% number of dots to skip
dot_skip = 2;
% polynomial order for z
order_z = 1;
% save calibration data? (true/false)
save_calibration_data = true;
% save camera model? (true/false)
save_camera_model = true;
% read previous calibration data? (true/false)
read_calibration_data = false;

% ==========================
%% run calibration
% ==========================
% load mask
mask = imread(fullfile('./sample-data/mask.tif'));
    
% create calibration directory
calibration_image_directory = fullfile('./sample-results/calibration/');
if ~exist(calibration_image_directory, 'dir')
    mkdir(calibration_image_directory);
end
    
% load sample reference image to be used for the calibration
im_ref = imread(fullfile('./sample-data/images/ordered/im_0002.tif'));
% mask image
im_ref_masked = create_masked_image(im_ref, mask);
% save masked image in the calibration folder
imwrite(im_ref_masked, fullfile(calibration_image_directory, 'im.tif'));

% ==========================
% load job
% ==========================
% load sample job file and adjust parameters
sample_job_filename = './sample-data/sample-calibration-job.mat';
sample_job = load(sample_job_filename);

% extract the calibration job
caljob = sample_job.datasave.caljob;

% call calibration code
calculate_camera_model_coefficients_reference_image_02(caljob, calibration_image_directory, save_calibration_data, save_camera_model, read_calibration_data, dot_diameter, dot_spacing, dot_skip, order_z);
