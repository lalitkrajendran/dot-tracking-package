function [Uc, Vc] = cross_correlate_dots_05(size_1, size_2, track, correlation_algorithm, subpixel_fit, zero_mean, min_sub)
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
%     rmin = min(r1,r2);
%     cmin = min(c1, c2);
%     
%     I1 = I1(1:rmin, 1:cmin);
%     I2 = I2(1:rmin, 1:cmin);
    
%     % ensure that the two windows are the same size by zero padding
%     I1 = padarray(I1, [rmax - r1, cmax - c1], 0, 'post');
%     I2 = padarray(I2, [rmax - r2, cmax - c2], 0, 'post');
    
%     r1_clip = floor(r1/2) * 2;
%     r2_clip = floor(r2/2) * 2;
%     c1_clip = floor(c1/2) * 2;
%     c2_clip = floor(c2/2) * 2;
%     
%     rmax = min(r1_clip, r2_clip);
%     cmax = min(c1_clip, c2_clip);
%     
%     I1 = I1(1:rmax, 1:cmax);
%     I2 = I2(1:rmax, 1:cmax);
    
    % location where the correlation estimate will be stored
%     X = round(0.5*size(I1,2));
%     Y = round(0.5*size(I1,1));
    X = round(0.5*min(size(I1,2), size(I2,2)));
    Y = round(0.5*min(size(I1,1), size(I2,1)));
    
    % -----------------------------------------------
    % processing settings for cross-correlation
    % -----------------------------------------------
    
    % window size
%     window_size_1 = size(I1);
%     window_size_2 = size(I2);
    window_size_1 = [size(I1,2), size(I1,1)];
    window_size_2 = [size(I2,2), size(I2,1)];
    
    % diameter for rpc filter
    filter_diameter = 2.8;

    % particle diameter
    particle_diameter = 0.5 * (size_1.XYDiameter(p_1,3) + size_2.XYDiameter(p_2,3));
    
%     % zero mean images?
%     if strcmp(correlation_algorithm, 'dcc')
%         zero_mean = 0;
%     else
%         zero_mean = 1;
%     end
%     zero_mean = 1;
    % subpixel estimator
    if strcmp(subpixel_fit, 'tpg')
        peak_locator = 1;
    elseif strcmp(subpixel_fit, 'fpg')
        peak_locator = 2;        
    elseif strcmp(subpixel_fit, 'lsg')
        peak_locator = 3;
    else
        fprintf('unknown subpixel fit. exiting.\n');
    end

    % find additional peaks?
    find_extrapeaks = 1;
    % apply a fractionally weighted filter?
    frac_filt = 0;
    % save cross correlation plane?
    saveplane = 1;
    % bulk velocities
    Ub = 0;
    Vb = 0;
    % zero pad images?
    zero_pad = 0;
    % effective window resolution
    window_resolution = [window_size_1; window_size_2];

    % -----------------------------------------------
    % cross-correlate image pair and estimate displacement from prana
    % -----------------------------------------------

%     [Xc,Yc,Uc,Vc,Cc,Dc,Cp, Dxc, Dyc] = PIVwindowed(I1,I2,correlation_algorithm,[window_size+2, window_size+2],window_resolution,zero_pad,filter_diameter*ones(1,2),zero_mean,peak_locator,find_extrapeaks,frac_filt,saveplane,X,Y,Ub,Vb);
%     [~,~,Uc,Vc,~,~,Cp,~,~] = PIVwindowed(I1,I2,correlation_algorithm,window_size_1,window_resolution,zero_pad,filter_diameter*ones(1,2),zero_mean,peak_locator,find_extrapeaks,frac_filt,saveplane,X,Y,Ub,Vb);
    
    if zero_mean
        I1 = I1 - mean(I1(:));
        I1(I1 < 0) = 0;
        I2 = I2 - mean(I2(:));
        I2(I2 < 0) = 0;
    end
    
    if min_sub
        I1 = I1 - min(I1(:));
        I2 = I2 - min(I2(:));
    end
    
    % calculate correlation coefficient
    Cp = conv2(I2, rot90(conj(I1),2));
%     Cp = xcorr2(I2, I1);
    
    % set negative values of the correlation plane to zero
%     Cp(Cp < 0) = 0;

%     % set negative values of the correlation plane to the inverse value
%     Cp(Cp < 0) = -Cp(Cp < 0);
    
    % offset correlation plane by its minimum value
    Cp = Cp - min(Cp(:));

    % subpixel estimation of displacement
%     [Uc,Vc,~,~]=subpixel_02(Cp,size(Cp,2),size(Cp,1),ones(size(Cp)),peak_locator,find_extrapeaks,[filter_diameter, filter_diameter], I1, I2);
%     Uc = Uc - size(I1,2);
%     Vc = Vc - size(I1,1);

    [Uc,Vc,~,~]=subpixel(Cp,size(Cp,2),size(Cp,1),ones(size(Cp)),peak_locator,find_extrapeaks,[particle_diameter, particle_diameter]);
    Uc = Uc - (-floor(size(Cp,2)/2)) + 1 - size(I1,2);
    Vc = Vc - (-floor(size(Cp,1)/2)) + 1 - size(I1,1);


%     if peak_locator == 1
%         [Uc,Vc,~,~,~,~]=Gaussfit(Cp, subpixel_fit, 4);
%     elseif peak_locator == 3
%         [Uc,Vc,~,~,~,~]=Leastsqrfit(Cp, subpixel_fit, 4);
%     end
    
    % if there are multiple peaks, then pick the one with displacement
    % closest to the tracking estimate
    if numel(Uc) > 1
        err_U = Uc - ((track(2) - c_min_2) - (track(1) - c_min_1));
        err_V = Vc - ((track(4) - r_min_2) - (track(3) - r_min_1));
        
        [~, minloc] = min(sqrt(err_U.^2 + err_V.^2));
        Uc = Uc(minloc);
        Vc = Vc(minloc);
    end    

%     [uc,vc,~,~,~,~]=Gaussfit(Cp,'tpg',4);
%     [uc,vc,~,~,~,~]=Leastsqrfit(Cp,'lsg',4);
%     uc = uc - size(I1,2);
%     vc = vc - size(I1,1);
    
%     % if one dimension of the correlation plane is biased, then account for
%     % a possible bias error in the subpixel fit
%     if mod(size(Cp,1),2) == 0
%         Vc = Vc - 1;
%     end
%     if mod(size(Cp,2),2) == 0
%         Uc = Uc - 1;
%     end
    
    
    % calculate final displacement in the global image co-ordinate system
    Uc_image = c_min_2 - c_min_1 + Uc;
    Vc_image = r_min_2 - r_min_1 + Vc;
%     Vc_image = r_min_1 - r_min_2 + Vc;
    
    % assign values to be returned by the function
    Uc = Uc_image;
    Vc = Vc_image;
%     Vc = -Vc_image;

end