function [p_matrix,peaks,num_p]=blob_segmentation(im,v)
%
%[p_matrix,peaks,num_p]=mapparticles(im,v)
%
%Simple particle identification method which, after global thresholding,
%groups clusters of adjacent pixels together ('blobs') and labels them as a
%particle.  Very fast and works well for sparse-to-medium image seeding
%densities.  The method has no mechanism to separate overlapped particles
%other tham aggressive thresholding, which adversly affects the image.
%
%INPUTS
%   im - image to be segmented (matrix-uint8)
%   v  - intial threshold value (num)
%
%OUTPUTS
%   p_matrix (2D array) - matrix of particle identification and extent
%   peaks (2D array) - matrix of identified image peaks (by image erosion)
%   num_p (num) - number of identified particles in p_matrix
%
%N.Cardwell - 2.8.2008 (modified from the work by M.Brady)
%N.Cardwell - 10.27.2009 (added a base thresholding operation and a peak
%   identification operation, removed section of code for the creation of
%   mapint and locxy - no longer needed with the TRACKING_V2 code)

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


%globally threshold the image
im_thresh=im;  im_thresh(im < v)=0;

%group adjacent pixels with nonzero intensity and assign each cluster a
%particle number and store in p_matrix, also returns total number of
%particles detected (num_p), searches by columns
[p_matrix,num_p]=bwlabel(im_thresh,8);

%determine the peak intensity pixel of each pixel blob
peaks=zeros(size(p_matrix));
for i=1:num_p
    [r_i,c_i]=find(p_matrix==i);  indx_i=sub2ind(size(peaks),r_i,c_i);
    I_i=im(indx_i);
    pixels_i=[r_i,c_i,I_i];
    pixels_i_sort=sortrows(pixels_i,3);
    peaks(pixels_i_sort(end,1),pixels_i_sort(end,2))=i;
end

end
