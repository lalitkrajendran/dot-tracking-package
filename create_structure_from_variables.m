function results_struct = create_structure_from_variables(varargin)
% Function to create a structure from the input variables
%
% INPUTS:
% varargin: any variable to be condensed into a structure
%
% OUTPUTS:
% results_struct: structure variable
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % calculate number of variables
    num_vars = numel(varargin);
    % declare results structure
    results_struct = struct;
    % loop through variables and assign values
    for var_index = 1:num_vars
        % extract variable name
        var_name = inputname(var_index);
        % extract variable value
        var_value = varargin{var_index};
        results_struct.(var_name) = var_value;
    end
end