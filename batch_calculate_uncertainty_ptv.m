clear
close all
clc

restoredefaultpath;

logical_string = {'False', 'True'};

% name for this processing case
case_name = 'id=apriori_size=lsg_corr=lsg';

current_image_directory = './sample-data/images/';
current_results_directory = './sample-results/';

% load processing parameters
params = load(fullfile(current_results_directory, case_name, 'id', 'id_params.mat'));
id = params.id;
params = load(fullfile(current_results_directory, case_name, 'size', 'size_params.mat'));
sizing = params.sizing;
params = load(fullfile(current_results_directory, case_name, 'tracks', 'track_params.mat'));
tracking = params.tracking;

% get directory listing for images
[image_files, num_image_files] = get_directory_listing(current_image_directory, 'ordered', 'im*.tif');
% get directory listing for size results
[size_files, num_size_files] = get_directory_listing(fullfile(current_results_directory, case_name, 'size'), 'file*.mat');
% get directory listing for tracks results
[track_files, num_track_files] = get_directory_listing(fullfile(current_results_directory, case_name, 'tracks'), 'file*.mat');        
% obtain list of files in the reference image directory
[ref_images, num_ref_images] = get_directory_listing(fullfile(current_image_directory, 'ref'), 'im*.tif');

% load mask
image_mask = load_image_mask('./sample-data/mask.tif');

% -------------------------------
%% extract id and sizing results for the tracked dots in the average image
% -------------------------------
% load first reference image
im_ref = imread(fullfile(image_files(2).folder, image_files(2).name));
im_ref = pre_process_image(im_ref, id.minimum_subtraction, id.minimum_intensity_level, image_mask);

% load sizing results for the first image pair
results = load(fullfile(size_files(1).folder, size_files(1).name));
size_ref = results.size_ref;

% ==========================================
%% identify dots on other reference images
% ==========================================
fprintf('dot id and sizing on reference images\n');
[ref_avg_props, ref_all_props] = calculate_average_reference_image_properties(im_ref, size_ref, ref_images, id, sizing, tracking, image_mask);
save(fullfile(current_results_directory, case_name, 'ref_images_stats.mat'), 'ref_avg_props', 'ref_all_props');

% =======================
% loop through all tracked files and calculate position uncertainty for the gradient image
% =======================
fprintf('calculating uncertainties\n');
for file_index = 1:num_track_files
    fprintf('file_index: %d\n', file_index);
    % ------------------
    % load processing results
    % ------------------
    % load tracks
    results = load(fullfile(track_files(file_index).folder, track_files(file_index).name));
    tracks = results.tracks;

    % calculate uncertainties
    uncertainty2D = calculate_uncertainty_dot_tracking_correlation(results.tracks, results.Cp, ref_avg_props, tracking);

    % append results to tracks
    save(fullfile(track_files(file_index).folder, track_files(file_index).name), 'uncertainty2D', '-append');
end
