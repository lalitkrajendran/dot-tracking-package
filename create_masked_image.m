function im2 = create_masked_image(im, mask)
% function to create a masked version of the image when the mask and the
% image are of different types. especially useful when the masked image is
% to be resaved in its original 
%
% INPUTS:
% im: image array
% mask: mask array
% 
% OUTPUTS:
% im2: masked image array with the same class as the original image
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
    
    % if image and the mask are the same integer class, multiply to create
    % the masked image
    if strcmp(class(im), class(mask))
        im2 = im .* mask;
    % if they are not the same class
    else
        % create a copy of the mask array
        mask2 = mask;
        % convert mask to a binary array (from 0s and 255s to 0s and 1s)
        mask2(mask2 > 0) = 1;
        % convert binary mask to the same integer class as the image to be
        % masked
        mask3 = cast(mask2, 'like', im);
        % multiply the image with the new binary mask of the same integer
        % class to create the masked image
        im2 = im .* mask3;
    end
end