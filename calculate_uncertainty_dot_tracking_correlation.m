function uncertainty2D = calculate_uncertainty_dot_tracking_correlation(tracks, Cp, ref_avg_props, tracking)
% Function to calculate position and displacement uncertainties from dot tracking
%
% INPUTS:
% tracks: tracked displacements
% Cp: correlation plane for each of the tracks
% ref_avg_props: average properties from the reference images
% tracking: tracking structure properties
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
    d_x = nans(num_tracks, 1);
    d_y = nans(num_tracks, 1);
    R = nans(num_tracks, 1);
    I = nans(num_tracks, 1);

    AR_U = nans(num_tracks, 1);
    AR_V = nans(num_tracks, 1);

    U_std = nans(num_tracks, 1);
    V_std = nans(num_tracks, 1);

    % =======================
    % loop through tracks and extract gradient dot properties
    % =======================
    parfor track_index = 1:num_tracks        
        % extract dot index 
        [p_ref, p_grad] = extract_dot_index(tracks(track_index, :));

        % calculate properties of the cross_correlation plane
        [~, D_l, I0_l, R, ~] = gauss_lsq_peakfit_2D_general_normalized(Cp{track_index}, tracking.correlation_correction.subpixel_fit, false);                            
        d_x(track_index) = D_l(1);
        d_y(track_index) = D_l(2);
        R(track_index) = R;
        I(track_index) = I0_l;        

        % ------------------
        % extract average properties
        % ------------------
        d_x_avg = ref_avg_props.d_x_ref_avg(p_ref);
        d_y_avg = ref_avg_props.d_y_ref_avg(p_ref);
        R_avg = ref_avg_props.R_ref_avg(p_ref);
        I_avg = ref_avg_props.I_ref_avg(p_ref);
        U_ref_std = ref_avg_props.U_ref_std(p_ref);
        V_ref_std = ref_avg_props.V_ref_std(p_ref);

        % ------------------
        % calculate amplification ratio                
        % ------------------
        [AR_U(track_index), AR_V(track_index)] = calculate_amplification_ratio(d_x_avg, d_y_avg, R_avg, I_avg, ...
                                                        d_x(track_index), d_y(track_index), R(track_index), I(track_index));

        % ------------------------------                
        %% calculate displacement uncertainty
        % ------------------------------                
        U_std(track_index) = AR_U(track_index) * U_ref_std;
        V_std(track_index) = AR_V(track_index) * V_ref_std;
    end
    
    % create results structure
    uncertainty2D = create_structure_from_variables(d_x, d_y, R, I, AR_U, AR_V, U_std, V_std);
end