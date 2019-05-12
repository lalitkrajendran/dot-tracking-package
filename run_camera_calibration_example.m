% example script to calculate mapping function coefficients from reference image

clear
close all
clc

restoredefaultpath;
addpath ../dot-tracking-package/

%% experiment settings
test_date = '2019-03-14';
% body placed in the tunnel
body = 'wedge';
% dot size on the target (mm)
dot_size = 0.15;
dot_size_string = num2str(['0_' num2str(dot_size*100) 'mm']);
% top level project directory
project_directory = '/scratch/shannon/c/aether/Projects/BOS/arms-supersonic-wind-tunnel/';
% top directory containing images for this case
top_image_directory = fullfile(project_directory, 'analysis', 'data', test_date, body, dot_size_string); 
% top directory to contain results for this case
top_results_directory = fullfile(project_directory, 'analysis', 'results', test_date, body, dot_size_string);

%% call calibration code
calibration_image_directory = fullfile(top_image_directory, 'calibration-crlb');
save_calibration_data = true;
save_camera_model = true;
read_calibration_data = false;
dot_diameter = dot_size * 1e-3;
dot_spacing = 2 * dot_size * 1e-3;
dot_skip = 3;
order_z = 1;

calculate_camera_model_coefficients_reference_image_02(calibration_image_directory, save_calibration_data, save_camera_model, read_calibration_data, dot_diameter, dot_spacing, dot_skip, order_z);