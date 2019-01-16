function [XYDiameter, mapsizeinfo, locxy, mapint]=combined_ID_size_apriori_02(im, x_ref, y_ref, d_p, sizing_method, min_area)
    % This function is modified from combined_partID.
    % It creates a result identical to the blob segmentation step 
    % but without any thresholding. It uses prior information about
    % location of dots on the pattern. 

    % INPUTS:
    % im - image
    % x_ref - x location of dots (pix.)
    % y_ref - y location of dots (pix.)
    % d_p - particle diameter (pix.)

    % modified by Lalit Rajendran. - 08/02/2018

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


    num_p = length(x_ref);

    p_matrix = zeros(size(im));
    p_ctr = 0;
    peaks=zeros(size(p_matrix));
    
    mapint = cell(1, num_p);
    locxy = zeros(num_p, 2);
    mapsizeinfo = zeros(num_p, 2);
    particleprops = nan * ones(num_p, 6);
    % create a particle intensity matrix
    for p = 1:num_p
        % array indices corresponding to the particle center
        r_p = floor(y_ref(p));
        c_p = floor(x_ref(p));

        % ignore if particle is not in field of view
        if r_p < 0 || c_p < 0 || r_p >= size(im,1) || c_p >= size(im,2)
            continue;
        end

        % extents of the particle image map
        r_min = floor(y_ref(p) - d_p/2);
        c_min = floor(x_ref(p) - d_p/2);

%         r_max = floor(y_ref(p) + d_p/2);
%         c_max = floor(x_ref(p) + d_p/2);
        r_max = ceil(y_ref(p) + d_p/2);
        c_max = ceil(x_ref(p) + d_p/2);

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

        if r_max - r_min <= 1 || c_max - c_min <= 1
            continue;
        end

        % intensity of the expected dot neighborhood
        im_p = im(r_min:r_max, c_min:c_max);
        mapint{p} = double(im_p);
        locxy(p,:) = [c_min, r_min];
        mapsizeinfo(p,:) = [r_max - r_min, c_max - c_min];
        
        %determine if blob_i has multiple peaks, if so then pass it to the
        %dynamic function for additional segmentation, otherwise keep the
        %particle and move on
        BW_max=imregionalmax(im_p);
        if nnz(BW_max)~=1
            %crop out the blob and segment it
            im_crop=im_p;
            [p_matrix_crop,peaks_crop,~]=dynamic_threshold_segmentation_v3(im_crop,0,0);

            % find number of elements for each peak
            num_pixels = zeros(1,max(peaks_crop(:)));
            stats = regionprops(p_mat
            for peak_index = 1:max(peaks_crop(:))
                num_pixels(peak_index) = sum(sum(p_matrix_crop == peak_index));
            end
            
            % find peak wiht maximum corresponding pixels
            [~, maxloc] = max(num_pixels);
            
            % set pixels with other values of peaks to zero
            im_p(p_matrix_crop ~= maxloc) = 0;
            mapint{p} = double(im_p);
            
%             [row, col, ~] = find(p_matrix_crop == maxloc);
%             locxy(p,:) = [c_min + min(col) - 1, r_min + min(row) - 1];
%             mapsizeinfo(p,:) = [max(row) - min(row), max(col) - min(col)];
        end

        %% sizing
        if sum(sum(im_p > 0)) < min_area
            particleprops(p,:)=[NaN,NaN,NaN,NaN,NaN,NaN];
            continue;
        end
        [x_cg,y_cg,D,I0,~,Meth]=Gaussfit(mapint{p},sizing_method,4);
        x_c = locxy(p,1)+x_cg-1;
        y_c = locxy(p,2)+y_cg-1;

        if nnz(isnan([x_c,y_c,D,I0])>1)
            if sizeprops.errors==0
                %retain the IWC estimate
            else
                particleprops(p,:)=[x_c,y_c,max(D),max(I0),p,Meth+1];
            end
        else
            particleprops(p,:)=[x_c,y_c,max(D),max(I0),p,Meth+1];
        end
    end
    
    XYDiameter = particleprops;
end

% 
% 
%     num_p = p_ctr;
% 
%     if strcmp(id_method, 'dynamic')
%         %preallocate arrays
%         p_matrix_new=zeros(size(p_matrix));  num_p_new=0;  peaks_new=zeros(size(p_matrix));
% 
%         %main program loop - to segment or not to segment, that is the question!
%         for p=1:num_p
%             blob_i=double(p_matrix==p);
% 
%             %use 'regionprops to det. the extent of blob_i and segment the blob out
%             STATS = regionprops(blob_i,'BoundingBox');
%             % account for empty bounding box
%             if isempty(STATS)
%                 continue;
%             end
%             B_Box=ceil(STATS.BoundingBox);
%             if B_Box(3)==1;  width=0;  else width=B_Box(3)-1;  end
%             if B_Box(4)==1;  height=0;  else height=B_Box(4)-1;  end
%             part_i=im( B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width );
%         %     part_i(part_i<v)=0;
% 
%             %determine if blob_i has multiple peaks, if so then pass it to the
%             %dynamic function for additional segmentation, otherwise keep the
%             %particle and move on
%             BW_max=imregionalmax(part_i);
%             if nnz(BW_max)~=1
%                 %crop out the blob and segment it
%                 im_crop=im(B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width);
%                 [p_matrix_crop,peaks_crop,num_p_crop]=dynamic_threshold_segmentation_v3(im_crop,0,0);
% 
%                 %reset the 'particle number' of the segmented blod so it can be
%                 %directly inserted into the 'new' parameters
%                 for j=1:num_p_crop
%                     num_p_new=num_p_new+1;
%                     peaks_crop(peaks_crop==j) = num_p_new;
%                     p_matrix_crop(p_matrix_crop==j) = num_p_new;
%                 end
% 
%                 %insert the segmented/cropped particle into the 'new' parameters
%                 peaks_new( B_Box(2):B_Box(2)+height ,...
%                     B_Box(1):B_Box(1)+width ) = peaks_crop;
%                 p_matrix_new( B_Box(2):B_Box(2)+height ,...
%                     B_Box(1):B_Box(1)+width ) = p_matrix_crop;
%             else
%                 %assign the particle to the 'new' parameters
%                 num_p_new=num_p_new+1;
%                 [r_peak,c_peak]=find(BW_max==1);
%                 peaks_new(B_Box(2)+(r_peak-1) , B_Box(1)+(c_peak-1)) = num_p_new;
%                 p_matrix_new(B_Box(2):B_Box(2)+height , B_Box(1):B_Box(1)+width) = num_p_new;
%             end
%         end
% 
%         %reassign output variables
%         num_p=num_p_new;
%         peaks=peaks_new;
%         p_matrix=p_matrix_new;
% 
%     end
%     
%     %% run dot sizing
%     particleprops = zeros(num_p, 6);
%     
%     if strcmpi(sizing_method,'TPG')%sizeprops.method==2
%         for p=1:num_p
%             [x_cg,y_cg,D,I0,~,Meth]=Gaussfit(mapint{p},sizing_method,4);
%             x_c = locxy(p,1)+x_cg-1;
%             y_c = locxy(p,2)+y_cg-1;
% 
%             if nnz(isnan([x_c,y_c,D,I0])>1)
%                 if sizeprops.errors==0
%                     %retain the IWC estimate
%                 else
%                     particleprops(p,:)=[x_c,y_c,max(D),max(I0),p,Meth+1];
%                 end
%             else
%                 particleprops(p,:)=[x_c,y_c,max(D),max(I0),p,Meth+1];
%             end
%         end
%     elseif strcmpi(sizing_method,'LSG') || strcmpi(sizing_method,'CLSG')%sizeprops.method==5 || sizeprops.method==6
%         for p=1:num_p
%             [x_cl,y_cl,D_l,I0_l,~,Meth] = Leastsqrfit(mapint{p},sizing_method,4);
%             x_c = locxy(p,1)+x_cl-1;
%             y_c = locxy(p,2)+y_cl-1;
% 
%             if nnz(isnan([x_c,y_c,D_l,I0_l]))>1
%                 if sizeprops.errors==0
%                     %retain the IWC estimate
%                 else
%                     particleprops(p,:)=[x_c,y_c,max(D_l),max(I0_l),p,Meth+2];
%                 end
%             else
%                 particleprops(p,:)=[x_c,y_c,max(D_l),max(I0_l),p,Meth+2];
%             end
%         end
%     end
%     
%     XYDiameter = particleprops;
%     
% end
