function [Uc, Vc] = cross_correlate_dots_07(im1, im2, size_1, size_2, track, subpixel_fit, zero_mean, min_sub)
% this function cross-correlates the intensity maps of two dots and
% estimates the displacement. The intensity maps are interrogation windows
% centered around a dot in the respective frames

% INPUTS:
% im1, im2 - image pair
% size_1, size_2 - sizing results
% track - tracking results
    
    % -----------------------------------------------
    % create interrogation windows in the two images
    % -----------------------------------------------
    
    % find dot numbers in the two frames
    p_1 = track(11);
    p_2 = track(12);
    
    % extract intensity maps
    I1 = size_1.mapint{p_1};
    I2 = size_2.mapint{p_2};
    
    % origin for intensity maps for each dot
    c_min_1 = size_1.locxy(p_1,1);
    r_min_1 = size_1.locxy(p_1,2);
    c_min_2 = size_2.locxy(p_2,1);
    r_min_2 = size_2.locxy(p_2,2);
    
    % find max size of the two windows
    [r1, c1] = size(I1);
    [r2, c2] = size(I2);    
    rmax = max(r1, r2);
    cmax = max(c1, c2);
        
    % extract windows from original image
    im1_window = im1(r_min_1:r_min_1+r1-1, c_min_1:c_min_1+c1-1);
    im2_window = im2(r_min_2:r_min_2+r2-1, c_min_2:c_min_2+c2-1);        
    
    % location where the correlation estimate will be stored
    X = round(0.5*min(size(I1,2), size(I2,2)));
    Y = round(0.5*min(size(I1,1), size(I2,1)));
    
    % -----------------------------------------------
    % processing settings for cross-correlation
    % -----------------------------------------------
    
    % window size
    window_size_1 = [size(I1,2), size(I1,1)];
    window_size_2 = [size(I2,2), size(I2,1)];
    
    % diameter for rpc filter
    filter_diameter = 2.8;

    % break if the sizing failed
    if isnan(size_1.XYDiameter(p_1, 3)) || isnan(size_2.XYDiameter(p_2, 3))
        Uc = NaN;
        Vc = NaN;
        return;
    else
        particle_diameter = 0.5 * (size_1.XYDiameter(p_1,3) + size_2.XYDiameter(p_2,3));
    end
    
    % assign subpixel estimator
    if strcmp(subpixel_fit, 'tpg')
        peak_locator = 1;
    elseif strcmp(subpixel_fit, 'fpg')
        peak_locator = 2;        
    elseif strcmp(subpixel_fit, 'lsg')
        peak_locator = 3;
    elseif strcmp(subpixel_fit, 'clsg')
        peak_locator = 4;
    else
        fprintf('unknown subpixel fit. exiting.\n');
    end

    % find additional peaks?
    find_extrapeaks = 1;

    % -----------------------------------------------
    % cross-correlate image pair and estimate displacement from prana
    % -----------------------------------------------
    % zero mean images if required
    if zero_mean
        I1 = I1 - mean(I1(:));
        I1(I1 < 0) = 0;
        I2 = I2 - mean(I2(:));
        I2(I2 < 0) = 0;
    end
    % perform minimum subtraction if required
    if min_sub
        I1 = I1 - min(im1_window(:)); %- min(I1(:));
        I1(I1 < 0) = 0;
        I2 = I2 - min(im2_window(:)); %- min(I2(:));
        I2(I2 < 0) = 0;
    end
    
    % calculate correlation coefficient
    Cp = conv2(I2, rot90(conj(I1),2));
        
    % offset correlation plane by its minimum value
    Cp = Cp - min(Cp(:));

    % subpixel estimation of displacement
    [Uc,Vc,~,~] = subpixel(Cp,size(Cp,2),size(Cp,1),ones(size(Cp)),peak_locator,find_extrapeaks,[particle_diameter, particle_diameter]);
    Uc = Uc - (-floor(size(Cp,2)/2)) + 1 - size(I1,2);
    Vc = Vc - (-floor(size(Cp,1)/2)) + 1 - size(I1,1);
    
    % if there are multiple peaks, then pick the one with displacement
    % closest to the tracking estimate
    if numel(Uc) > 1
        err_U = Uc - ((track(2) - c_min_2) - (track(1) - c_min_1));
        err_V = Vc - ((track(4) - r_min_2) - (track(3) - r_min_1));
        
        [~, minloc] = min(sqrt(err_U.^2 + err_V.^2));
        Uc = Uc(minloc);
        Vc = Vc(minloc);
    end       
    
    % calculate final displacement in the global image co-ordinate system
    Uc_image = c_min_2 - c_min_1 + Uc;
    Vc_image = r_min_2 - r_min_1 + Vc;
    
    % assign values to be returned by the function
    Uc = Uc_image;
    Vc = Vc_image;

end