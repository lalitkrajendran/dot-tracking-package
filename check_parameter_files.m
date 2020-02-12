function [io, id, sizing, tracking] = check_parameter_files(io, id, sizing, tracking)
% function to ensure that the parameter structure has all the required
% fields
    if ~isfield(tracking, 'validation')
        tracking.validation = false;
    end

    if ~isfield(id, 'minimum_subtraction')
        id.minimum_subtraction = false;
    end
end