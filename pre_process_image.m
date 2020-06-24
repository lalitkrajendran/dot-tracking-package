function im_new = pre_process_image(im, perform_minimum_subtraction, minimum_intensity_level,image_mask)
% This function applies common pre-processing operations for the image
% 1. Converts to double
% 2. Peforms minimum subtraction if required
% 3. Flips image upside down
%
% INPUTS:
% im: image
% id: data structure containing 
    % convert images to double arrays
    im_new = double(im);
    
    % peform minimum subtraction if required
    if perform_minimum_subtraction
        im_new = im_new - minimum_intensity_level;
        im_new(im_new < 0) = 0;
    end
    
    % flip images upside down
    im_new = flipud(im_new);
    
    % image masking
    im_new = im_new .* image_mask;
end