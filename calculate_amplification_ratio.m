function [AR_x, AR_y] = calculate_amplification_ratio(d_x_ref, d_y_ref, R_ref, I_ref, d_x_grad, d_y_grad, R_grad, I_grad)
% This function calculates the uncertainty amplification ratio for a dot in
% a BOS image as described in the following paper.
%
% Rajendran et. al. (2020). Uncertainty amplification due to 
% density/refractive- index gradients in volumetric PTV and BOS experiments
%
% INPUT:
% _ref: dot properties in the reference image
% _grad: dot properties in the gradient image
% d_x, d_y: dot diameter along x and y [pix.]
% R: correlation coefficient 
% I: peak intensity 
%
% OUTPUT:
% AR_x, AR_y: uncertainty amplification ratios along x and y
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
    
    % calculate amplification ratio along x
    AR_x = sqrt(d_x_grad/d_x_ref) * sqrt(d_y_ref/d_y_grad) *...
        I_ref/I_grad * ...
        ((1 - R_grad^2)/(1 - R_ref^2))^0.25;
    
    % calculate amplification ratio along y
    AR_y = sqrt(d_y_grad/d_y_ref) * sqrt(d_x_ref/d_x_grad) *...
        I_ref/I_grad * ...
        ((1 - R_grad^2)/(1 - R_ref^2))^0.25;

end