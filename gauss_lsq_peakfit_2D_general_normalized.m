function [x_c, d_p, I0, R, covariance_matrix] = gauss_lsq_peakfit_2D_general_normalized(I, sizing_method, display_figure)

% This function does Gaussian least square fitting of 2D particle image
% function and evaluates the full covariance matrix to resolve any rotation
% in corr-elation function
% Input:
% I: particle image
% sizing_method: 'lsg', 'clsg'
% display_figure: option to display intermediate results of the
% optimization
% Output:
% x_c: particle centroid
% d_p: particle diameter
% I0: peak intensity
% R: correlation coefficient
% covariance_matrix: The estimated Covariance for Gaussian least square fit

% This function is written by Sayantan Bhattacharya

% ----------------------------------
% Improvements to previous versions
% ----------------------------------

% use a normalized Gaussian model. See the following paper for more
% information.

% Rajendran et. al. (2020) Uncertainty amplification due to
% density/refractive-index gradients in Background-Oriented Schlieren
% experiments.

% it uses a normalized version of the precision matrix such that it is
% directly the inverse of the covariance matrix

% it returns the x and y diameters instead of major and minor axes

% allows the centroids to also be variable during the
% optimization process

% allows both 'lsg' and 'clsg' to be used

% Modified by Lalit Rajendran (lrajendr@purdue.edu)

% ----------------------------------


    %Subtracting the minimum of the correlation plane
    I = I - min(I(:));

    % Max Value of correlation plane
    Imax = max(I(:));

    % This is the linear index of the correlation maximum (this only takes the first maximum
    % value if multiple are found - this is unlikely, but should be addressed)
    % max_linear_index=find(I==max(I(:)),1);
    [ii_max, jj_max] = find(I==Imax, 1);

    % calculate centroid and diaemter using intensity weighted centroid
    [x_centroid,y_centroid,diameter,~] = centroidfit(I, [0, 0]);
    centroid = [x_centroid, y_centroid];
    diameter = [diameter, diameter];

    % extract centroids
    iwc_flag = 0;
    
    % check if any of the results are nan
    if ~isnan(centroid(1))
        jj_max = centroid(1);
    else
        iwc_flag = iwc_flag + 1;    
    end

    if ~isnan(centroid(2))
        ii_max = centroid(2);
    else
        iwc_flag = iwc_flag + 1;    
    end
    
    % account for case when centroidfit failed in diameter calculation but could calculate a
    % centroid
    if isnan(diameter(1)) && iwc_flag == 0
        diameter(1) = size(I,2);
    end

    if isnan(diameter(2)) && iwc_flag == 0
        diameter(2) = size(I,1);
    end
    
    % if diameter is greater than the size of the image, then fix that problem
    if diameter(1) > size(I, 2)
        diameter(1) = size(I, 2);
    end

    if diameter(2) > size(I, 1)
        diameter(2) = size(I, 1);
    end

    % exit function if iwc fails
    if iwc_flag > 0
        fprintf('IWC fit failed. Exiting.\n');
        x_c = ones(1, 2) * NaN;
        d_p = ones(1, 2) * NaN;
        I0 = NaN;
        theta = NaN;
        CovarianceMat = ones(2) * NaN;
        return;
    end

    % extract diameters
    Dx = diameter(1);% X correlation diameter for jj direction
    Dy = diameter(2);% Y correlation diameter for ii direction

    % assign window to be used for least squares fitting
    x_min = 1;
    x_max = size(I, 2);
    y_min = 1;
    y_max = size(I, 1);

    % Extracted peak region
    peak_region=double(abs(I(y_min:y_max,x_min:x_max)));

    % Fitting options
    options=optimset('MaxIter', 20,'MaxFunEvals', 50,'TolX', 1e-9,'TolFun', 1e-9,...
        'Display','off','DiffMinChange',1e-7,'DiffMaxChange',1e-3,...
        'Algorithm','Levenberg-Marquardt', 'Scaleproblem', 'jacobian'); %, 'FinDiffType', 'central');
    % options = optimoptions(@lsqnonlin,'MaxIterations',1200,'MaxFunctionEvaluations',5000,'Display','off');


    %Grid points (This is double checked and matches the cropped region for
    %least square fit
    [xg, yg]=meshgrid(x_min:x_max,y_min:y_max);

    locxy_i=[yg(:),xg(:)];
    Xp = locxy_i(:,2);
    Yp = locxy_i(:,1);

    Xvec=[Xp(:)';Yp(:)'];

    % initial conditions: x=[I0, Clevel, xc, yc, a, c, b]
    x0=double([Imax, 0, x_centroid, y_centroid, Dx, Dy, 0]);
    % lower bound (L-M optimization cannot handle LB/UB)
    LB = [];
    % upper bound
    UB = [];

    % Run solver; default to iwc if it fails
    try
        % run optimization
        [xvars, ~, ~, exitflag]=lsqnonlin(@calculate_intensity_profile, x0, LB, UB, options, ...
                                            peak_region, Xvec, sizing_method, display_figure);

        % extract optimization results if it was succesful
        if exitflag >= 0 %== 1 || exitflag == 3
            % ----------------------------
            % extract optimization results
            % ----------------------------
            % peak intensity
            I0 = xvars(1);
            % x centroid
            xloc = xvars(3); %jj_max;
            % y centroid
            yloc = xvars(4); %ii_max;
            % x diameter
            Dx = xvars(5); %xvars(3);
            % y diameter
            Dy = xvars(6); %xvars(4);
            % correlation coefficient
            R = xvars(7);
            % if x diameter is larger than size of the image, fix this
            if Dx > size(I, 2)
                Dx = size(I, 2);
            end

            % if y diameter is larger than size of the image, fix this
            if Dy > size(I, 1)
                Dy = size(I, 1);
            end

            % construct the covariance matrix
            covariance_matrix =[(Dx/4)^2, R * Dx/4 * Dy/4;...
                            R * Dx/4 * Dy/4, (Dy/4)^2];            
            % flag to indicate optimization was successful
            fit_flag = 0;
        else
            % flag to indicate optimization was unsuccessful
            fit_flag = 1;
        end

    catch
        % flag to indicate optimization was unsuccessful        
        fit_flag = 1;
    end
    
    % use IWC results if case of failed optimization
    if fit_flag == 1
        % inform the user
        fprintf('\n Least squared fit failed, defaulting to IWC \n');
        
        % x centroid
        xloc=centroid(1);
        % y centroid
        yloc=centroid(2);
        % x diameter
        Dx=diameter(1);
        % y diameter
        Dy=diameter(2);
        % peak intensity
        I0=Imax;
        % correlation coefficient
        R = 0;
        % covariance matrix
        covariance_matrix =[(Dx/4)^2, R * Dx/4 * Dy/4;...
                        R * Dx/4 * Dy/4, (Dy/4)^2];
    end
    
    % aggreate centroids
    x_c = [xloc, yloc];
    % aggregate diameters
    d_p = [Dx, Dy];
                
end


function display_intensity_profile(I_meas, I_fit)
% This function displays the measured intensity profile from the image and
% an intensity profile from a Gaussian approximation

    % extract max measured intensity
    I_max = max(I_meas(:));
    
    figure(100)
    % display measured intensity
    subplot(2, 3, 1)
    imagesc(I_meas)
    colormap(flipud(gray))
    caxis([0 I_max])
    title('Data')

    % display predicted intensity
    subplot(2, 3, 2)
    imagesc(I_fit)
    colormap(flipud(gray))
    caxis([0 I_max])
    title('Fit')

    % display error
    subplot(2, 3, 3)
    imagesc(abs(I_fit - I_meas))
    colormap(flipud(gray))
    caxis([0 I_max])
    title('Error')
    
    % display measured intensity, calculated intensity and error along a
    % horizontal slice through the center
    subplot(2, 3, 4)
    r0 = round(size(I_meas, 1)/2);
    plot(I_meas(r0,:), '*')
    hold on
    plot(I_fit(r0,:), '-')
    plot(abs(I_fit(r0,:) - I_meas(r0,:)), 'o')
    hold off
    title(['r = ' num2str(r0, '%d')])
    
    % display measured intensity, calculated intensity and error along a
    % vertical slice through the center
    subplot(2, 3, 5)
    c0 = round(size(I_meas, 2)/2);
    plot(I_meas(:, c0), '*')
    hold on
    plot(I_fit(:, c0), '-')
    plot(abs(I_fit(:, c0) - I_meas(:, c0)), 'o')
    hold off
    title(['c = ' num2str(c0, '%d')])
    l1 = legend('Data', 'Fit', 'Error', 'location', 'eastoutside');
    l1.Position = [0.7425    0.2367    0.0926    0.1099];
    set(gcf, 'Position', [122   252   993   619]);
    
    % pause figure for viewer
    pause(0.1);
end

function F = calculate_intensity_profile(x, mapint_i, Xvec, sizing_method, display_figure)
% This function calculates an intensity profile for the Gaussian given the
% dot parameters

    % -------------------------------------
    % extract variables from input argument
    % -------------------------------------
    % peak intensity
    I0 = x(1);
    % background intensity level
    Clevel = x(2);
    % x centroid
    x_c = x(3); 
    % y centroid
    y_c = x(4); 
    % x standard deviation (1/4th of the diameter)
    eta_x = x(5)/4;
    % y standard deviation (1/4th of the diameter)
    eta_y = x(6)/4;
    % correlation coefficient
    R = x(7);
    
    % construct the covariance matrix
    C = [eta_x^2, R * eta_x * eta_y; R * eta_x * eta_y, eta_y^2];
    
    % calculate the precision matrix
    B = inv(C);
    
    % assumed function for the intensity profile
    fun = @(x, y) exp(-0.5 * ((x - x_c) .* (B(1,1) * (x - x_c) + B(1,2) * (y - y_c)) + ...
                              (y - y_c) .* (B(2,1) * (x - x_c) + B(2,2) * (y - y_c))  ...
                              ));                          

    % initialize intensity profile of the dot
    gauss_int = zeros(size(Xvec,2),1);
    
    % loop through pixels and calculate intensity profile based on current
    % guesses for the dot parameters
    for i=1:size(Xvec,2)
        % current pixel
        xi = Xvec(1,i);
        yi = Xvec(2,i);
        
        % use the formulation for the discrete least squares
        if strcmpi(sizing_method, 'lsg')
            gauss_int(i) = fun(xi, yi);
        % use the formulation for the continuous least squares            
        elseif strcmpi(sizing_method, 'clsg')
            gauss_int(i) = integral2(fun, xi-0.5, xi+0.5, yi-0.5, yi+0.5);
        end
    end
    
    % account for background intensity and scaling
    gauss_int = Clevel + I0 * gauss_int;
    
    % calculate error in the guessed intensity profile
    F = mapint_i(:) - gauss_int;
        
    % display intensity maps if required
    if display_figure
        display_intensity_profile(mapint_i, reshape(gauss_int, size(mapint_i, 1), size(mapint_i, 2)));
    end    
end
