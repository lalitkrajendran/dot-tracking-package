function XYDiameter = perform_dot_sizing(num_p, mapint, locxy, sizing_method, default_iwc)
    
    particleprops = nan * ones(num_p, 8);        
    
    for p = 1:num_p
        %% sizing
        if nnz(isnan(locxy(p, :))) > 0
            continue;
        end
        % estimate centroid using IWC either if it is specified or if a
        % default estimate is to be computed        
        if strcmp(sizing_method, 'iwc') || default_iwc
            [particleprops(p,1),particleprops(p,2),particleprops(p,3),particleprops(p,5)]=...
                 centroidfit(mapint{p},locxy(p,:));
                particleprops(p,1) = particleprops(p,1) - 1;
                particleprops(p,2) = particleprops(p,2) - 1;
                % assign y diameter to be same as x diameter
                particleprops(p,4) = particleprops(p,3);
                % peak intensity
                particleprops(p,6) = particleprops(p,5);
                % correlation coefficient
                particleprops(p,5) = 0;
                % particle index
                particleprops(p,7) = p;
                % sizing method index
                particleprops(p,8) = 1;
        end
        
        % -----------------------------------------------------------------
        % estimate centroid using a Gaussian fit
        % -----------------------------------------------------------------
        if strcmp(sizing_method, 'tpg')
            % 3 PT Gaussian fit
            [x_cg,y_cg,D,I0,~,Meth]=Gaussfit(mapint{p},sizing_method,4);
            x_c = locxy(p,1)+x_cg-1;
            y_c = locxy(p,2)+y_cg-1;

            if nnz(isnan([x_c,y_c,D(1),D(2),max(I0)])) < 1
                particleprops(p,:)=[x_c, y_c, D(1), D(2), 0, max(I0), p, Meth+1];
            end
            
        elseif strcmpi(sizing_method,'LSG') || strcmpi(sizing_method,'CLSG')
            % Least Squares Gaussian fit
            [x_cl, D_l, I0_l, R, ~] = gauss_lsq_peakfit_2D_general_normalized(mapint{p}, sizing_method, false);
            % calculate centroid location on global co-ordinate system
            x_c = locxy(p,1) + x_cl(1) - 1;
            y_c = locxy(p,2) + x_cl(2) - 1;

            % Taking care of peak rotation
            Dx = D_l(1); %sqrt( (cos(alpha)^2 * D_l(1)^2 + sin(alpha)^2 * D_l(2)^2) );
            Dy = D_l(2); %sqrt( (sin(alpha)^2 * D_l(1)^2 + cos(alpha)^2 * D_l(2)^2) );
            
            if nnz(isnan([x_c,y_c,max(D_l),max(I0_l)])) < 1
                particleprops(p,:)=[x_c,y_c,Dx,Dy,R,I0_l,p,6];
            else
                fprintf('LSG fit failed\n');
            end

        end        
    end
    
    XYDiameter = particleprops;
end