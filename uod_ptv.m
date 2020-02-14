function [u_val, v_val, val, r0_u, r0_v] = uod_ptv(x, y, u, v, replace_vectors, residual_threshold, epsilon_uod)
% Function to perform universal outlier detection for unstructured PTV data
% based on the method by Duncan et. al.
% Duncan, J., Dabiri, D., Hove, J., & Gharib, M. (2010). Universal outlier detection for particle image velocimetry (PIV) and particle tracking velocimetry (PTV) data. Measurement Science and Technology, 21(5), 057002.
%
% INPUTS:
% x,y: co-ordiantes of particles in the first frame
% u,v: velocities of the corresponding particles
% replace_vectors: option to replace invalid vectors by interpolation of
%                  surrounding vectors (true/false)
% residual_threshold: threshold to classify a vector as invalid. Duncan et.
% al. recommend a range of 2-4 with a smaller threshold retaining less
% spurious vectors (stricter)
% epsilon_uod: threshold for the minimum allowable displacement based on
% Westerweel and Scarano (2005). Recommended: 0.1 pix.
%
% OUTPUTS:
% u_val, v_val: validated velocities
% val: array signifying whether a vector was valid of invalid (1 for
% invalid, 0 for valid)
% r0_u, r0_v: normalized residuals for each vector
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    if nargin < 5
        replace_vectors = false;
    elseif nargin < 6
        residual_threshold = 2;
        epsilon_uod = 0.1;
    elseif nargin < 7
        epsilon_uod = 0.1;
    end

    % identify neighbors using Delaunay triangulation
    num_tracks = size(x,1);

    % starting co-ordinates of tracks
    P = [x, y];
    
    % find neighborhood of tracks using delaunay triangulation
    DT = delaunayTriangulation(P);

    % find all attached edges that are the neighboring particles to the
    % given particle
    E = edges(DT);

    u_val = ones(size(u)) * NaN;
    v_val = ones(size(v)) * NaN;
    val = zeros(size(v));
    r0_u = ones(size(u)) * NaN;
    r0_v = ones(size(v)) * NaN;
    
    for track_index = 1:num_tracks

        %% identify neighbors

        % find the IDs of the points that make up the edge with the current
        % point
        neighbor_ID = [E(E(:,1) == track_index,2); E(E(:,2) == track_index, 1)];
        
        if isempty(neighbor_ID)
            continue;
        end
        % extract co-ordinates of the current point
        x0 = DT.Points(track_index,1);
        y0 = DT.Points(track_index,2);

        % extract co-ordinates of the neighboring points
        xi = DT.Points(neighbor_ID,1);
        yi = DT.Points(neighbor_ID,2);

        % calculate distance of each neighbor point to the current point
        di = sqrt((xi - x0).^2 + (yi - y0).^2);

        % extract velocity of current point
        U0 = u(track_index);
        V0 = v(track_index);
        
        % extract velocity of neighboring points
        Ui = u(neighbor_ID);
        Vi = v(neighbor_ID);
        
        %% calculate tolerance
        
        % calculate median of all distances to neighboring points
        b = median(di);
        
        % calculate discriminant of quadratic equation
        D = sqrt(b^2 + 4 * epsilon_uod);
        
        % calculate tolerance
        epsilon_a = 0.5 * (-b + D);

        %% calculate normalized residual
        
        % ---------------
        % U
        % ---------------
        % numerator
        numerator = abs(U0/(b + epsilon_a) - median(Ui./(di + epsilon_a)));
        % denominator
        denominator = median(Ui./(di + epsilon_a) - median(Ui./(di + epsilon_a))) + epsilon_a;
        % residual
        r0_u(track_index) = numerator/denominator;

        % ---------------
        % V
        % ---------------
        % numerator
        numerator = abs(V0/(b + epsilon_a) - median(Vi./(di + epsilon_a)));
        % denominator
        denominator = median(Vi./(di + epsilon_a) - median(Vi./(di + epsilon_a))) + epsilon_a;
        % residual
        r0_v(track_index) = numerator/denominator;
        
        %% remove residual and replace by interpolated value if required 
        
        % if residual is greater than threshold identify track to be a bad vector
        % and replace by interpolation of neighbors

        if r0_u(track_index) > residual_threshold
            val(track_index) = 1;
            if replace_vectors
                F_U = scatteredInterpolant(xi, yi, Ui, 'natural');
                interp_val = F_U(x0, y0);
                % if the interpolated value is not empty, then assign it to
                % the track
                if ~isempty(interp_val) && ~isnan(interp_val)
                    u_val(track_index) = interp_val;
                else
                    u_val(track_index) = NaN;
                end
            else
                u_val(track_index) = NaN;
            end
        else
            u_val(track_index) = u(track_index);
        end
        
        if r0_v(track_index) > residual_threshold
            val(track_index) = 1;
            if replace_vectors
                F_V = scatteredInterpolant(xi, yi, Vi, 'natural');
                % if the interpolated value is not empty, then assign it to
                % the track
                interp_val = F_V(x0, y0);
                if ~isempty(interp_val) && ~isnan(interp_val)
                    v_val(track_index) = interp_val;
                else
                    v_val(track_index) = NaN;
                end
            else
                v_val(track_index) = NaN;                
            end
        else
            v_val(track_index) = v(track_index);
        end
        
    end    

end
