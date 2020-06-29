function ref_avg_props = calculate_average_reference_image_properties(im_ref, size_ref, ref_images, id, sizing, tracking, image_mask)
 % Function to calculate average properties from a series of reference images
 %
 % INPUTS:
 % im_ref: reference image used in tracking
 % size_ref: sizing results for the reference image used in tracking
 % ref_images: listing of all reference images
 % id: dot identification settings
 % sizing: dot sizing settings
 % tracking: dot tracking settings
 % image_mask: image mask
 %
 % OUTPUTS:
 % ref_avg_props: data structure containing average properties of the result
 % of correlating each reference image with the one used for tracking.
 %
 % AUTHOR:
 % Lalit Rajendran (lrajendr@purdue.edu)   
 
    % number of reference images
    num_ref_images = numel(ref_images);

    % number of dots
    num_p = size(size_ref.XYDiameter, 1);

    % declare cells and arrays to hold co-ordinates of identified dots
    size_ref_all = cell(1, num_ref_images);
    U_ref_all = nans(num_p, num_ref_images);
    V_ref_all = nans(num_p, num_ref_images);
    d_x_ref_all = nans(num_p, num_ref_images);
    d_y_ref_all = nans(num_p, num_ref_images);
    R_ref_all = nans(num_p, num_ref_images);
    I_ref_all = nans(num_p, num_ref_images);

    % loop through all images
    for image_index = 1:num_ref_images    
        fprintf('image: %d\n', image_index);
        
        % load image
        im = imread(fullfile(ref_images(image_index).folder, ref_images(image_index).name));    

        % prepare image for processing
        im = pre_process_image(im, id.minimum_subtraction, id.minimum_intensity_level, image_mask);
        
        % copy locxy and mapsize info
        size_current.locxy = size_ref.locxy;
        size_current.mapsizeinfo = size_ref.mapsizeinfo;

        % extract dot intensity maps
        size_current.mapint = extract_dot_intensity_map(im, size_current.locxy, size_current.mapsizeinfo);

        % dot sizing
        size_current.XYDiameter = perform_dot_sizing(num_p, size_current.mapint, size_current.locxy, sizing.centroid_subpixel_fit, sizing.default_iwc);

        % create track array
        tracks_ref = zeros(num_p, 12);
        tracks_ref(:, 11) = 1:num_p;
        tracks_ref(:, 12) = 1:num_p;
        
        % cross correlate dots and size correlation plane
        parfor p = 1:num_p
            % cross correlate dots
            [U_ref_all(p, image_index), V_ref_all(p, image_index), Cp] = cross_correlate_dots_07(im_ref, im, size_ref, size_current, tracks_ref(p, :), ...
                        tracking.correlation_correction.subpixel_fit, tracking.correlation_correction.zero_mean, tracking.correlation_correction.min_sub);                        
            
            % calculate properties of the cross_correlation plane
            [~, D_l, I0_l, R, ~] = gauss_lsq_peakfit_2D_general_normalized(Cp, tracking.correlation_correction.subpixel_fit, false);                            
            d_x_ref_all(p, image_index) = D_l(1);
            d_y_ref_all(p, image_index) = D_l(2);
            R_ref_all(p, image_index) = R;
            I_ref_all(p, image_index) = I0_l;
        end
    end

    % calculate average properties
    d_x_ref_avg = mean(d_x_ref_all, 2, 'omitnan');
    d_y_ref_avg = mean(d_y_ref_all, 2, 'omitnan');
    R_ref_avg = mean(R_ref_all, 2, 'omitnan');
    I_ref_avg = mean(I_ref_all, 2, 'omitnan');
    
    % calculate standard deviation of reference image position
    U_ref_std = std(U_ref_all, [], 2, 'omitnan');
    V_ref_std = std(V_ref_all, [], 2, 'omitnan');

    % create structure with these properties
    ref_avg_props = create_structure_from_variables(d_x_ref_avg, d_y_ref_avg, R_ref_avg, I_ref_avg, U_ref_std, V_ref_std);
end