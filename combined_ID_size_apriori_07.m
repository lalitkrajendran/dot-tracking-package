function [XYDiameter, peaks, mapsizeinfo, locxy, mapint]=combined_ID_size_apriori_07(im, x_ref, y_ref, d_p, sizing_method, default_iwc, min_area, W_area, W_intensity, W_distance)
    % This function was modifed from combined_partID and particle_sizing
    % written by Cardwell.
    % This function performs particle/dot identification and sizing without
    % the use of intensity thresholds for the segmentation process. It 
    % instead uses prior information about location of dots on the image.

    % INPUTS:
    % im - image
    % x_ref - x location of dots [pix.]
    % y_ref - y location of dots [pix.]
    % d_p - particle diameter [pix.]
    % sizing_method - subpixel fit (e.g., iwc/tpg/lsg)
    % default_iwc - option to default to the iwc estimate if the gaussian fits fail.
    % min_area - minimum area for a blob to be considered a particle [pix.]
    % W_area - area weight for calculating the peak coefficient used to
    % differentiate between true and false blobs
    % W_intensity - weight for peak intensity
    % W_distance - weight for distance from center
    
    % OUTPUTS:
    %   XYDiameter = 6-column matrix; 1st column is x-centorid location, 2nd
    %       column is y-centroid location, 3rd column is particle diameter, 4th
    %       column is true max. intensity, 5th column is particle id#, 6th
    %       column is sizing method
    %   mapsizeinfo - (2 column array) defines the size (row col) of the each
    %       sized particle
    %   locxy - (2 column array) ROW/COL? location associated with the upper 
    %       left pixel of the particle's rectangular projection

    % modified by Lalit Rajendran. - 09/13/2018

    %This function uses a combination of the 'blob' and 'dynamic threshold
    %segmentation' algorithms in an attempt to reduce the computational cost of
    %the dynamic while still retaining its ability to accuartely segment 
    %overlaped particles.  Special care must be taken in the intital
    %thresholding of the image to partially segment the blobs without removing
    %too many of the low intensity particle images
    %
    %(beta) N.Cardwell - 10.28.2009

    %     This file is part of prana, an open-source GUI-driven program for
    %     calculating velocity fields using PIV or PTV.
    %     Copyright (C) 2012  Virginia Polytechnic Institute and State
    %     University
    % 
    %     prana is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    % 
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    % 
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.


    % number of dots expected to be in the image
    num_p = length(x_ref);

    % ensure that the diameter is an array to avoid dimension mismatch
    if numel(d_p) == 1
        d_p = d_p * ones(num_p, 1);
    end
    
    % intialize arrays
    p_matrix = zeros(size(im));
    p_ctr = 0;
    peaks=zeros(size(p_matrix));    
    mapint = cell(1, num_p);
    locxy = nan * ones(num_p, 2);
    mapsizeinfo = nan * ones(num_p, 2);
    particleprops = nan * ones(num_p, 6);
    % create a particle intensity matrix
    for p = 1:num_p
        %% 1. create window around expected dot locatoin
        
        % image array indices corresponding to the particle center
        r_p = floor(y_ref(p));
        c_p = floor(x_ref(p));

        if isnan(r_p) || isnan(c_p)
            continue;
        end
        
        % ignore if particle is not in field of view
        if r_p <= 0 || c_p <= 0 || r_p >= size(im,1) || c_p >= size(im,2)
            continue;
        end

        % extents of the particle image map
        r_min = floor(y_ref(p) - d_p(p)/2);
        c_min = floor(x_ref(p) - d_p(p)/2);

        r_max = ceil(y_ref(p) + d_p(p)/2);
        c_max = ceil(x_ref(p) + d_p(p)/2);

        % clip boundaries
        if r_min < 1
            r_min = 1;
        end
        if c_min < 1
            c_min = 1;
        end
        if r_max > size(im,1)
            r_max = size(im,1);
        end
        if c_max > size(im,2)
            c_max = size(im,2);
        end

        % ignore if resulting image map is too small
        if r_max - r_min <= 1 || c_max - c_min <= 1
            continue;
        end

        % intensity of the expected dot neighborhood
        im_p = im(r_min:r_max, c_min:c_max);
        mapint{p} = double(im_p);
        locxy(p,:) = [c_min, r_min];
        mapsizeinfo(p,:) = [r_max - r_min + 1, c_max - c_min + 1];
        
        % if all pixels are zero, then continue
        if(sum(im_p) == 0)
            continue;
        end
        
        %% 2. Perform dynamic segmentation to detect additional peaks
        % crop out the blob and segment it
        im_crop=im_p;
%         [p_matrix_crop,peaks_crop,~]=dynamic_threshold_segmentation_v3(im_crop,0,0);
        [p_matrix_crop,peaks_crop,~]=dynamic_threshold_segmentation_v5(im_crop,0,0, r_p, c_p);
        flag = 1;
        loop_ctr = 0;
        % check if multiple peaks have been detected
        if max(peaks_crop(:) > 1)            
            %% combine peaks that are close to each other
            while flag > 0 && loop_ctr <=20
                flag = 0;
                loop_ctr = loop_ctr + 1;
                
                % find peak locations
                [r_peak,c_peak] = ind2sub(size(im_p), find(peaks_crop > 0));            

                % ensure no duplicate peaks
                peak_indices = peaks_crop(peaks_crop(:) > 0);            
                if length(unique(peak_indices)) < length(peak_indices)
                    [~, ind] = unique(peaks_crop(peaks_crop>0));
                    duplicate_ind = setdiff(1:size(peaks_crop(peaks_crop>0), 1), ind);
                    peaks_crop(r_peak(duplicate_ind), c_peak(duplicate_ind)) = 0;

                    % find new peak locations if duplicates were identified and
                    % removed
                    [r_peak,c_peak] = ind2sub(size(im_p), find(peaks_crop > 0));
                    peak_indices = unique(peak_indices);
                end

                % find distance between peaks
                distance_matrix = ones(length(r_peak)) * NaN;
                truth_matrix = zeros(length(r_peak));
                for i = 1:length(r_peak)
                    for j = i+1:length(r_peak)
                        distance_matrix(i, j) = sqrt((r_peak(i)-r_peak(j)).^2 + (c_peak(i) - c_peak(j)).^2);
                        distance_matrix(j, i) = distance_matrix(i,j);                    
                    end

                    % find nearest neighbor
                    [min_distance(i), neighbor(i)] = min(distance_matrix(i,:));    
                end

                % pair peaks if distance is less than dot radius
                paired = zeros(1,length(r_peak));
                for i = 1:length(r_peak)
                    % check if distance between peaks is less than dot radius
                    if ~paired(i) && ~paired(neighbor(i)) && min_distance(i) < d_p(p)/2
                        flag = flag + 1;
                        % pair peaks
                        paired(i) = 1;
                        paired(neighbor(i)) = 1;

                        % delete original peaks
                        peaks_crop(r_peak(i), c_peak(i)) = 0;
                        peaks_crop(r_peak(neighbor(i)), c_peak(neighbor(i))) = 0;

                        % add new peak at midpoint to be same as ith peak
                        mid_point = [0.5*(r_peak(i) + r_peak(neighbor(i))), 0.5*(c_peak(i) + c_peak(neighbor(i)))];
                        peaks_crop(round(mid_point(1)), round(mid_point(2))) = peak_indices(i);

                        % combine p_matrix for the two peaks
                        p_matrix_crop(p_matrix_crop == peak_indices(neighbor(i))) = peak_indices(i);
                    end                
                end

                % re-number peaks
                peak_nums_sorted = sort(unique(peaks_crop(peaks_crop > 0)));
                for peak_index = 1:length(peak_nums_sorted)
                    peaks_crop(peaks_crop == peak_nums_sorted(peak_index)) = peak_index;
                    p_matrix_crop(p_matrix_crop == peak_nums_sorted(peak_index)) = peak_index;
                end

            end
        
            %% detect the peak that corresponds to the particle
            
            % number of peaks identified from dynamic segmentation
            num_peaks = max(peaks_crop(:));

            % arrays to hold peak properties
            peak_area = zeros(1, num_peaks);
            peak_intensities = zeros(1, num_peaks);
            peak_distance = zeros(1, num_peaks);

            % calculate statistics of the identified blobs
            stats = regionprops(p_matrix_crop, double(im_p), 'area', 'maxintensity', 'centroid');
            if length(stats) < max(peaks_crop(:))
                keyboard;
            end
            for peak_index = 1:max(peaks_crop(:))
                % find number of elements for each peak
                peak_area(peak_index) = stats(peak_index).Area; %sum(sum(p_matrix_crop == peak_index));
                if peak_area(peak_index) == 0
                    continue;
                end
                % find peak intensity
                peak_intensities(peak_index) = stats(peak_index).MaxIntensity; %im_p(peaks_crop == peak_index);               
                % calculate peak distance from predicted position
                peak_distance(peak_index) = (stats(peak_index).Centroid(1) - (x_ref(p) - c_min + 0.5)).^2 + (stats(peak_index).Centroid(2) - (y_ref(p) - r_min + 0.5)).^2;
%                 peak_distance(peak_index) = (stats(peak_index).Centroid(1) - (size(im_crop,2)/2+0.5)).^2 + (stats(peak_index).Centroid(2) - (size(im_crop,1)/2+0.5)).^2;
            end

            % ignore this dot is all identified blobs have areas that are
            % too small
            if max(peak_area) < min_area
                continue;
            end

            % calculate a weighted coefficient of the peaks
            peak_coefficients = (W_area * peak_area/max(peak_area) + ...
                W_intensity * peak_intensities/max(peak_intensities) + ...
                W_distance * (1 - peak_distance/max(peak_distance))) * 1/(W_area + W_intensity + W_distance);

            % sort peaks in descending order of the weighted coefficient
            [~, i_coeff] = sort(peak_coefficients, 'descend');
            % pick the peak with the highest coefficient to be the correct
            % blob
            maxloc = i_coeff(1);

            % set pixels with other values of peaks to zero
            peak_loc = find(peaks_crop == maxloc);
            [r_peak, c_peak] = ind2sub(size(p_matrix_crop), peak_loc);

            % fill in the peaks
            peaks(r_min + r_peak-1, c_min + c_peak-1) = p;
            
            % set pixels corresponding to other peaks as 0
            im_p(p_matrix_crop ~= maxloc) = 0;
            
            % ensure that the dot image is convex (to avoid incursions of
            % zero pixels in the dot intensity map in the dynamic
            % segmentation process)
            
            % first create binary image
            im_bw = double(im_p);
            im_bw(im_bw > 0) = 1;
            
            % identify the convex hull in the image and the bounding box
            stats = regionprops(im_bw, double(im_p), 'conveximage', 'boundingbox');
            
            % identify extents of the bounding box
            rmin_crop = ceil(stats.BoundingBox(2));
            cmin_crop = ceil(stats.BoundingBox(1));
            rmax_crop = ceil(stats.BoundingBox(2)) + stats.BoundingBox(4) - 1;
            cmax_crop = ceil(stats.BoundingBox(1)) + stats.BoundingBox(3) - 1;
            
            % ensure that the image is convex and without holes
            im_p = double(im_crop(rmin_crop:rmax_crop, cmin_crop:cmax_crop)) .* stats.ConvexImage;
            mapint{p} = double(im_p);
            
            % update the dot window location on the global image
            r_min = r_min + rmin_crop - 1;
            c_min = c_min + cmin_crop - 1;
                        
            locxy(p,:) = [c_min, r_min];
            mapsizeinfo(p,:) = [stats.BoundingBox(4), stats.BoundingBox(3)];

%             % crop im to only contain pixels corresponding to the actual
%             % dot
%             [r,c] = ind2sub(size(p_matrix_crop), find(p_matrix_crop == maxloc));
%             im_p_2 = im_p(min(r):max(r), min(c):max(c));
%             im_p = im_p_2;
%             mapint{p} = double(im_p);
%             
%             % update origin location on the global co-ordinate system
%             if min(r) > 1 || min(c) > 1
%                 r_min = r_min + min(r) - 1;
%                 c_min = c_min + min(c) - 1;
%                 r_max = r_min + max(r) - 1;
%                 c_max = c_min + max(c) - 1;
% 
%                 locxy(p,:) = [c_min, r_min];
%                 mapsizeinfo(p,:) = [r_max - r_min, c_max - c_min];
%             end
            
        else
            peaks(r_p, c_p) = p;
        end

        %% sizing

        % estimate centroid using IWC either if it is specified or if a
        % default estimate is to be computed        
        if strcmp(sizing_method, 'iwc') || default_iwc
            [particleprops(p,1),particleprops(p,2),particleprops(p,3),particleprops(p,4)]=...
                 centroidfit(mapint{p},locxy(p,:));
                particleprops(p,1)=particleprops(p,1) - 1;
                particleprops(p,2)=particleprops(p,2) - 1;
                particleprops(p,5)=p;
                particleprops(p,6)=1;
        end
        
        % -----------------------------------------------------------------
        % estimate centroid using a Gaussian fit
        % -----------------------------------------------------------------
        if strcmp(sizing_method, 'tpg')
            % 3 PT Gaussian fit
            [x_cg,y_cg,D,I0,~,Meth]=Gaussfit(mapint{p},sizing_method,4);
            x_c = locxy(p,1)+x_cg-1;
            y_c = locxy(p,2)+y_cg-1;

%             if nnz(isnan([x_c,y_c,D,I0])) < 1
%                 particleprops(p,:)=[x_c,y_c,max(D),max(I0),p,Meth+1];
%             end
            if nnz(isnan([x_c,y_c,max(D),max(I0)])) < 1
                particleprops(p,:)=[x_c,y_c,max(D),max(I0),p,Meth+1];
            end
            
        elseif strcmpi(sizing_method,'LSG') || strcmpi(sizing_method,'CLSG')
            % Least Squares Gaussian fit
            [x_cl,y_cl,D_l,I0_l,~,Meth] = Leastsqrfit(mapint{p},sizing_method,4);
            x_c = locxy(p,1)+x_cl-1;
            y_c = locxy(p,2)+y_cl-1;

%             if nnz(isnan([x_c,y_c,D_l,I0_l])) < 1
%                 particleprops(p,:)=[x_c,y_c,max(D_l),max(I0_l),p,Meth+2];
%             end
            if nnz(isnan([x_c,y_c,max(D_l),max(I0_l)])) < 1
                particleprops(p,:)=[x_c,y_c,max(D_l),max(I0_l),p,Meth+2];
            else
                fprintf('LSG fit failed\n');
            end

        end        
    end
    
    XYDiameter = particleprops;
end
