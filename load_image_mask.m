function image_mask = load_image_mask(image_mask_filepath)
    % load image mask
    image_mask = imread(image_mask_filepath);
    % convert mask to 0 and 1
    image_mask(image_mask > 0) = 1;
    % convert mask to double
    image_mask = double(image_mask);
    % flip mask upside down
    image_mask = flipud(image_mask);
end