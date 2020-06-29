function mapint = extract_dot_intensity_map(im, locxy, mapsizeinfo)        
% Function to extract the intensity map of a dot
%
% INPUTS:
% im: image
% locxy: upper left corner of each dot sub image
% mapsizeinfo: size of each dot sub image
%
% OUTPUTS:
% mapint: intensity map for each dot
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)    
    
    % number of dots
    num_p = size(locxy, 1);
    % initialize array to hold intensity maps
    mapint = cell(1, num_p);
    % loop through all the dots
    for p = 1:num_p
        % extract extents of the dot
        rmin = locxy(p, 2);
        cmin = locxy(p, 1);
        rmax = rmin + mapsizeinfo(p, 1) - 1;
        cmax = cmin + mapsizeinfo(p, 2) - 1;

        % if the dot had failed identification, then skip
        if nnz(isnan([rmin, cmin, rmax, cmax])) < 1
            mapint{p} = im(rmin:rmax, cmin:cmax);
        else
            mapint{p} = nans(mapsizeinfo(1), mapsizeinfo(2));
        end
    end
end