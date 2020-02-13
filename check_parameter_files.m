function [io, id, sizing, tracking] = check_parameter_files(io, id, sizing, tracking)
% function to ensure that the parameter structure has all the required
% fields
    %% check i/o parameters
    if ~isfield(io, 'display_intermediate_progress')
        io.display_intermediate_progress = false;
    end

    %% check id parameters
    if ~isfield(id, 'minimum_subtraction')
        id.minimum_subtraction = false;
    end

    %% check tracking parameters
    if ~isfield(tracking, 'validation')
        tracking.validation = false;
    end
    
    if ~isfield(tracking.validation, 'uod_epsilon')
        tracking.validation.uod_epsilon = 0.1;
    end
end