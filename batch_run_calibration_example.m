% example script to calculate mapping function coefficients from reference image

clear
close all
clc

restoredefaultpath;
addpath ../dot-tracking-package/
addpath /scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/

% ==========================
%% calibration settings
% ==========================
% dot size (mm)
dot_size = 0.042;
dot_diameter = dot_size * 1e-3;
dot_spacing = 2 * dot_size * 1e-3;
dot_skip = 2;
order_z = 1;
save_calibration_data = true;
save_camera_model = true;
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
im_ref = imread(fullfile('./sample-data/images/im_0002.tif'));
% mask image
im_ref_masked = create_masked_image(im_ref, mask);
% save masked image in the calibration folder
imwrite(im_ref_masked, fullfile(calibration_image_directory, 'im.tif'));

% call calibration code
calculate_camera_model_coefficients_reference_image_02(calibration_image_directory, save_calibration_data, save_camera_model, read_calibration_data, dot_diameter, dot_spacing, dot_skip, order_z);
