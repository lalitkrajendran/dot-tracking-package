function size_results_new = remove_nan_entries_from_sizing(size_results)
% Function to remove nan entries from sizing results
%
% INPUTS:
% size_results: results from sizing
% 
% OUTPUTS:
% size_results_new: updated results with nan entries removed
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
    
    % identify rows containing nan elements
    [nan_r, ~] = find(isnan(size_results.XYDiameter));
    % create data structure to hold new results
    size_results_new = size_results;
    % remove corresponding entries
    size_results_new.XYDiameter(nan_r, :) = [];
    size_results_new.mapsizeinfo(nan_r, :) = [];
    size_results_new.locxy(nan_r, :) = [];
    size_results_new.mapint(nan_r) = [];
end