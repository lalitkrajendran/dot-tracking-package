% example script to calculate mapping function coefficients from reference image

clear
close all
clc

restoredefaultpath;
addpath ../dot-tracking-package/
addpath /scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/

%% experiment settings

% date on which the test was conducted
date = '2019-03-22';
% electrode gap (mm)
electrode_gap_array = [2, 5];
% imaging views
views = {'front-view'; 'side-view'; 'top-view'};

%% read/write settings

% top directory containing images
top_image_directory = '/scratch/shannon/c/aether/Projects/plasma-induced-flow/analysis/data/dbd/';
% top directory containing results
top_results_directory = '/scratch/shannon/c/aether/Projects/plasma-induced-flow/analysis/results/dbd/';

%% calibration settings
% dot size (mm)
dot_size = 0.042;
dot_diameter = dot_size * 1e-3;
dot_spacing = 2 * dot_size * 1e-3;
dot_skip = 2;
order_z = 1;
save_calibration_data = true;
save_camera_model = true;
read_calibration_data = false;
%% run calibration

% =====================================
% loop through electrode gaps
% =====================================
for electrode_gap = electrode_gap_array
    fprintf('Gap: %d mm\n', electrode_gap);
    current_gap_directory = fullfile(top_image_directory, date, [num2str(electrode_gap) 'mm']);
    % =====================================
    % loop through views
    % =====================================
    for view_index = 1:numel(views)
        fprintf('%s\n', views{view_index});
        current_view_directory = fullfile(current_gap_directory, views{view_index});
        % get list of folders containing different pulse parameters
        pulse_parameter_folders = dir(fullfile(current_view_directory, 'bos', '*V*ns'));
        % load mask
        mask = imread(fullfile(current_view_directory, 'bos', 'mask.tif'));
        % =====================================
        % loop through pulse parameters
        % =====================================
        for pulse_parameter_index = 1:length(pulse_parameter_folders)
            fprintf('case: %s\n', pulse_parameter_folders(pulse_parameter_index).name);
            current_pulse_parameter_directory = fullfile(current_view_directory, 'bos', pulse_parameter_folders(pulse_parameter_index).name);
            
            % create calibration directory
            calibration_image_directory = fullfile(current_pulse_parameter_directory, 'calibration');
            if ~exist(calibration_image_directory, 'dir')
                mkdir(calibration_image_directory);
            end
            
            % load sample reference image to be used for the calibration
            im_ref = imread(fullfile(current_pulse_parameter_directory, 'shot01', 'ref', 'im__0001.tif'));
            % mask image
            im_ref_masked = create_masked_image(im_ref, mask);
            % save masked image in the calibration folder
            imwrite(im_ref_masked, fullfile(calibration_image_directory, 'im.tif'));
            
            % call calibration code
            calculate_camera_model_coefficients_reference_image_02(calibration_image_directory, save_calibration_data, save_camera_model, read_calibration_data, dot_diameter, dot_spacing, dot_skip, order_z);
        end
    end
end
