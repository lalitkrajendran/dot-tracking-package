function uncertainty2D = calculate_uncertainty_dot_tracking(tracks, size_ref, size_grad, ref_avg_props)
% Function to calculate position and displacement uncertainties from dot tracking
%
% INPUTS:
% tracks: tracked displacements
% size_ref, size_grad: id and sizing results from ref and grad images
% ref_avg_props: average properties from the reference images
%
% OUTPUTS:
% uncertainty2D: structure containing the position and displacement uncertainties
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)    

    % calculate number of tracks
    num_tracks = size(tracks, 1);
    
    % ------------------
    % initialize variables
    % ------------------
    X_ref_std_tracked = nans(num_tracks, 1);
    Y_ref_std_tracked = nans(num_tracks, 1);

    X_grad_std_tracked = nans(num_tracks, 1);
    Y_grad_std_tracked = nans(num_tracks, 1);

    AR_x_ref = nans(num_tracks, 1);
    AR_y_ref = nans(num_tracks, 1);

    AR_x_grad = nans(num_tracks, 1);
    AR_y_grad = nans(num_tracks, 1);

    U_std = nans(num_tracks, 1);
    V_std = nans(num_tracks, 1);

    % =======================
    % loop through tracks and extract gradient dot properties
    % =======================
    parfor track_index = 1:num_tracks
        % extract dot index 
        [p_ref, p_grad] = extract_dot_index(tracks(track_index, :));
        
        % ------------------
        % extract dot properties
        % ------------------
        % reference
        ref_props = extract_dot_properties(size_ref.XYDiameter, p_ref);
        % gradient
        grad_props = extract_dot_properties(size_grad.XYDiameter, p_grad);

        % ------------------
        % extract average properties
        % ------------------
        d_x_avg = ref_avg_props.d_x_ref_avg(p_ref);
        d_y_avg = ref_avg_props.d_y_ref_avg(p_ref);
        R_avg = ref_avg_props.R_ref_avg(p_ref);
        I_avg = ref_avg_props.I_ref_avg(p_ref);
        X_ref_std = ref_avg_props.X_ref_std(p_ref);
        Y_ref_std = ref_avg_props.Y_ref_std(p_ref);

        % ------------------
        % calculate amplification ratio                
        % ------------------
        % reference
        [AR_x_ref(track_index), AR_y_ref(track_index)] = calculate_amplification_ratio(d_x_avg, d_y_avg, R_avg, I_avg, ...
                                                                                        ref_props.d_x, ref_props.d_y, ref_props.R, ref_props.I);
        % gradient
        [AR_x_grad(track_index), AR_y_grad(track_index)] = calculate_amplification_ratio(d_x_avg, d_y_avg, R_avg, I_avg, ...
                                                                                        grad_props.d_x, grad_props.d_y, grad_props.R, grad_props.I);

        % ------------------------------                
        %% calculate position uncertainty
        % ------------------------------                
        % reference
        X_ref_std_tracked(track_index) = X_ref_std * AR_x_ref(track_index);
        Y_ref_std_tracked(track_index) = Y_ref_std * AR_y_ref(track_index);
        
        % gradient
        X_grad_std_tracked(track_index) = X_ref_std * AR_x_grad(track_index);
        Y_grad_std_tracked(track_index) = Y_ref_std * AR_y_grad(track_index);                                                                                                        

        % ----------------------------------
        % calculate displacement uncertainty
        % ----------------------------------
        U_std(track_index) = sqrt(X_ref_std_tracked(track_index).^2 + X_grad_std_tracked(track_index).^2);
        V_std(track_index) = sqrt(Y_ref_std_tracked(track_index).^2 + Y_grad_std_tracked(track_index).^2);
    end
    
    % create results structure
    uncertainty2D = create_structure_from_variables(AR_x_ref, AR_y_ref, AR_x_grad, AR_y_grad, ...
                                                    X_ref_std_tracked, Y_ref_std_tracked, X_grad_std_tracked, Y_grad_std_tracked, ...
                                                    U_std, V_std);
end