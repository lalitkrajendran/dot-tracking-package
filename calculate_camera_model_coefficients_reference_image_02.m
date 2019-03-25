%%%%% DRAFT CODE, NOT COMPLETE %%%

clear
close all
clc

dbstop if error

read_calibration_data = false;

% top_image_directory = fullfile('/home/shannon/c/aether/Projects/BOS/spark-induced-flow/analysis/data/2018-08-20/', '90_450_0_1');
% top_image_directory = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-05/0_25mm/allpsi/';
% top_image_directory = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-11/0_15mm/allpsi/';
% top_image_directory = '/home/shannon/c/aether/Projects/BOS/nozzle/analysis/data/images/2018-10-11/0_1mm/allpsi/';
% top_image_directory = '/home/shannon/c/aether/Projects/BOS/arms-supersonic-wind-tunnel/analysis/data/2019-02-25/wedge-front/';
top_image_directory = '/scratch/shannon/c/aether/Projects/BOS/arms-supersonic-wind-tunnel/analysis/data/2019-03-14/wedge/0_25mm/';

save_calibration_data = true;
save_camera_model = true;
save_filepath = fullfile(top_image_directory, 'calibration');

if read_calibration_data
    save_filename = 'calibration_data.mat';
    temp_data = load(fullfile(save_filepath, save_filename));
    calibration_data = temp_data.calibration_data;
    calibration_plane_data = temp_data.calibration_plane_data;
    calibration_data.y_pixel_number = 309;
else
    % load sample job file and adjust parameters
    sample_job_filename = '/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/sample-job-files/sample-stereo-calibration-job.mat';
    sample_job = load(sample_job_filename);

    % extract the calibration job
    caljob = sample_job.datasave.caljob;

    % set calibration images
%     calibration_image_directory = fullfile(top_image_directory, 'test1', 'processing', 'average-ref-cropped');
%     calibration_image_directory = fullfile(top_image_directory, 'test1', 'processing', 'single-ref');
    calibration_image_directory = fullfile(top_image_directory, 'calibration');
    [files, ~] = get_directory_listing(calibration_image_directory, 'im*.tif');

    num_files = 1;
    num_cameras = 1;

    caljob.JOBFILE.Camera_Numbers_List = 1:num_cameras;

    caljob.calimagelist = cell(2, num_files);
    for t=1:num_files
        caljob.calimagelist{1,t} = fullfile(calibration_image_directory, files(t).name);
        caljob.calimagelist{2,t} = fullfile(calibration_image_directory, files(t).name);
    end

    caljob.JOBFILE.CalImageList = caljob.calimagelist;
    caljob.JOBFILE.Plane_Numbers_List = 1:num_files;
    caljob.JOBFILE.Z_Grid_Start = 0;
    caljob.JOBFILE.Plane_Spacing = 0;
    caljob.JOBFILE.Grid_Point_Diameter = 1; % 20e-3;
    caljob.JOBFILE.X_Grid_Spacing = 2*caljob.JOBFILE.Grid_Point_Diameter; %2*180e-3;
    caljob.JOBFILE.Y_Grid_Spacing = 2*caljob.JOBFILE.Grid_Point_Diameter; %2*180e-3;

    % load a sample image
    im = imread(fullfile(files(1).folder, files(1).name));
    [NR, NC] = size(im);
    caljob.x_pixel_number = NC;
    caljob.y_pixel_number = NR;
    
    % select control points
    [calibration_data,calibration_plane_data]=camera_calibration_new(caljob.JOBFILE);

    if save_calibration_data
        % save the object and image co-ordinates of the points on the calibration
        % grid
        save_filename = 'calibration_data.mat';
        save(fullfile(save_filepath, save_filename), 'calibration_data', 'calibration_plane_data');
    end
end

%%
% extract co-ordinates
allx1data(:,1)   = calibration_data.x_world_full{1};        % contains all x,y,z data for camera 1
allx1data(:,2)   = calibration_data.y_world_full{1};
allx1data(:,3)   = calibration_data.z_world_full{1};

allX1data(:,1)   = calibration_data.x_image_full{1};        % contains all X,Y data for camera 1
allX1data(:,2)   = calibration_data.y_image_full{1};

% create properties for a a fake second camera to run the fitmodels
% function
allx2data(:,1)   = calibration_data.x_world_full{1};        % contains all x,y,z data for camera 2
allx2data(:,2)   = calibration_data.y_world_full{1};
allx2data(:,3)   = calibration_data.z_world_full{1};    

allX2data(:,1)   = calibration_data.x_image_full{1};
allX2data(:,2)   = calibration_data.y_image_full{1};

% number of pixels
rA1=NR;

% flip y coordinate
% plus 1 is because the first index is at 1, not 0
% allX1data(:,2)  = rA1-allX1data(:,2) + 1; 
% allX2data(:,2)  = rA1-allX2data(:,2) + 1;

%%
% specify camera model
order_z   = 1; % cubic xy, quadratic z
optionsls   = [];

% fit camera model
[a_cam1, a_cam2, aXcam1, aYcam1, aXcam2, aYcam2, convergemessage] = fitmodels(allx1data,...
    allx2data,allX1data,allX2data,order_z,optionsls);
if save_camera_model
    save_filename = ['camera_model_type=' num2str(order_z) '.mat'];
    save(fullfile(save_filepath, save_filename), 'a_cam1', 'a_cam2', 'aXcam1', 'aYcam1', 'aXcam2', 'aYcam2');
end

