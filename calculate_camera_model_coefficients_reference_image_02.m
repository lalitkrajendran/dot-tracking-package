function calculate_camera_model_coefficients_reference_image_02(calibration_image_directory, save_calibration_data, save_camera_model, read_calibration_data, dot_diameter, dot_spacing, dot_skip, order_z) 
% function to calculate the camera model for a bos experiment from the
% reference image
%
% INPUTS:
% calibration_image_directory: directory containing calibration images
% save_calibration_data: save detected grid points? (true/false)
% save_camera_model: save camera model coefficients? (true/false)
% read_calibration_data: read previously stored calibration data?
% (true/false)
% dot_diameter: diameter of a dot on the bos target (m)
% dot_spacing: spacing between adjacent on the target (m)
% dot_skip: number of dots you are going to skip in the identification
% (default = 0)
% order_z: order of the z polynomial in the mapping function coefficient.
% (1 or 2)

% OUTPUS:
% None. The results will be saved in a mat file in the calibration
% directory.

    save_filepath = calibration_image_directory;

    if read_calibration_data
        save_filename = 'calibration_data.mat';
        temp_data = load(fullfile(save_filepath, save_filename));
        calibration_data = temp_data.calibration_data;
        calibration_plane_data = temp_data.calibration_plane_data;
        calibration_data.y_pixel_number = 309;
    else
        % ==========================
        % load job
        % ==========================
        % load sample job file and adjust parameters
        sample_job_filename = '/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/sample-job-files/sample-stereo-calibration-job.mat';
        sample_job = load(sample_job_filename);

        % extract the calibration job
        caljob = sample_job.datasave.caljob;

        % set calibration images
        [files, ~] = get_directory_listing(calibration_image_directory, 'im*.tif');

        % ==========================
        % set up calibration data structure
        % ==========================
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
        caljob.JOBFILE.Grid_Point_Diameter = dot_diameter * 1e3; %0.25; % 20e-3;
        caljob.JOBFILE.X_Grid_Spacing = (1 + dot_skip) * dot_spacing * 1e3; %2*4*caljob.JOBFILE.Grid_Point_Diameter; %2*180e-3;
        caljob.JOBFILE.Y_Grid_Spacing = (1 + dot_skip) * dot_spacing * 1e3; %2*4*caljob.JOBFILE.Grid_Point_Diameter; %2*180e-3;

        % load a sample image
        im = imread(fullfile(files(1).folder, files(1).name));
        [NR, NC] = size(im);
        caljob.x_pixel_number = NC;
        caljob.y_pixel_number = NR;

        % select control points
        [calibration_data, calibration_plane_data] = camera_calibration_new(caljob.JOBFILE);

        if save_calibration_data
            % save the object and image co-ordinates of the points on the calibration
            % grid
            save_filename = 'calibration_data.mat';
            save(fullfile(save_filepath, save_filename), 'calibration_data', 'calibration_plane_data');
        end
    end

    % ==========================
    % extract co-ordinates
    % ==========================
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
    
    % specify camera model
    optionsls = [];

    % fit camera model
    [a_cam1, a_cam2, aXcam1, aYcam1, aXcam2, aYcam2, convergemessage] = fitmodels(allx1data, allx2data, allX1data, allX2data, ...
                                                                                    order_z, optionsls);
    if save_camera_model
        save_filename = ['camera_model_type=' num2str(order_z) '.mat'];
        save(fullfile(save_filepath, save_filename), 'a_cam1', 'a_cam2', 'aXcam1', 'aYcam1', 'aXcam2', 'aYcam2');
    end

end
