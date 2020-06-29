function [p_ref, p_grad] = extract_dot_index(current_track)
% Function to extract dot index from tracking results
%
% INPUTS:
% current_track: row from the tracks array
%
% OUTPUTS:
% p_ref: dot index in reference image
% p_grad: dot index in gradient image
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
    
    % reference
    p_ref = current_track(12);        
    % gradient
    p_grad = current_track(11);
end