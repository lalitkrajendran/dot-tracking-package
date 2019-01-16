function [xp, yp] = calculate_predicted_dot_positions_02(X1, Y1, initialization_method, X, Y, U, V)
% this function calculates the predicted positions of dots in frame 2 based
% on their position in frame 1 and the displacement field obtained from
% cross-correlation
    
    if strcmp(initialization_method, 'rays')
        % create a scattered interpolant
        F_U = scatteredInterpolant(X, Y, U);
        F_V = scatteredInterpolant(X, Y, V);
        u_guess = F_U(X1, Y1);
        v_guess = F_V(X1, Y1);
    elseif strcmp(initialization_method, 'correlation')
        % interpolate to find the guess for tracking
        u_guess = interp2(X, Y, U, X1, Y1, 'spline');
        v_guess = interp2(X, Y, V, X1, Y1, 'spline');
    end
%     if size(X,2) > 1
%         % interpolate to find the guess for tracking
%         u_guess = interp2(X, Y, U, X1, Y1, 'spline');
%         v_guess = interp2(X, Y, V, X1, Y1, 'spline');
%     else
%         % interpolate to find the guess for tracking
%         u_guess = interp1(X, Y, U, X1, Y1, 'spline');
%         v_guess = interp1(X, Y, V, X1, Y1, 'spline');
%     end
    
    % remove nan elements
    u_guess(isnan(u_guess)) = 0;
    v_guess(isnan(v_guess)) = 0;

    % calculate predicted positions
    xp = X1 + u_guess;
    yp = Y1 + v_guess;
end